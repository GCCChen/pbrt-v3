
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// core/film.cpp*
#include "film.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"
#include "sbf/VectorNf.h"
#include "sbf/SBFCommon.h"
#include "sbf/TwoDArray.h"
#include "sbf/CrossBilateralFilter.h"
#include <iostream>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);

// Film Method Definitions
Film::Film(const Point2i &resolution, const Bounds2f &cropWindow,
           std::unique_ptr<Filter> filt, Float diagonal,
           const std::string &filename, Float scale,
         const vector<float> &interParams,
         const vector<float> &finalParams,
         float sigmaN, float sigmaR, float sigmaD,
         float interMseSigma, float finalMseSigma,
         Float maxSampleLuminance): fullResolution(resolution),
      diagonal(diagonal * .001),
      filter(std::move(filt)),
      rFilter(filter),
      interParams(interParams),
      finalParams(finalParams),
      sigmaN(sigmaN),
      sigmaR(sigmaR),
      sigmaD(sigmaD),
      interMseSigma(interMseSigma),
      finalMseSigma(finalMseSigma),
      filename(filename),
      scale(scale),
      maxSampleLuminance(maxSampleLuminance) {
    // Compute film image bounds
    croppedPixelBounds =
        Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                         std::ceil(fullResolution.y * cropWindow.pMin.y)),
                 Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                         std::ceil(fullResolution.y * cropWindow.pMax.y)));
    LOG(INFO) << "Created film with full resolution " << resolution <<
        ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
        croppedPixelBounds;

    // Allocate film image storage
    pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
    filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

    int xPixelCount = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
    int yPixelCount = croppedPixelBounds.pMax.y - croppedPixelBounds.pMin.y;

    colImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    varImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    featureImg = TwoDArray<Feature>(xPixelCount, yPixelCount);
    featureVarImg = TwoDArray<Feature>(xPixelCount, yPixelCount);

    norImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    rhoImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    depthImg = TwoDArray<float>(xPixelCount, yPixelCount);
    rhoVarImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    norVarImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    depthVarImg = TwoDArray<float>(xPixelCount, yPixelCount);
    
    fltImg = TwoDArray<Color>(xPixelCount, yPixelCount);
    minMseImg = TwoDArray<float>(xPixelCount, yPixelCount);
    adaptImg = TwoDArray<float>(xPixelCount, yPixelCount);
    sigmaImg = TwoDArray<Color>(xPixelCount, yPixelCount);

}

Bounds2i Film::GetSampleBounds() const {
    Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
                               Vector2f(0.5f, 0.5f) - filter->radius),
                         Ceil(Point2f(croppedPixelBounds.pMax) -
                              Vector2f(0.5f, 0.5f) + filter->radius));
    return (Bounds2i)floatBounds;
}

Bounds2f Film::GetPhysicalExtent() const {
    Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
    return std::unique_ptr<FilmTile>(new FilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void Film::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int i = 0; i < 3; i++) {
            pixel.Lrgb[i] = pixel.sqLrgb[i] =
                pixel.normal[i] = pixel.sqNormal[i] = 
                pixel.rho[i] = pixel.sqRho[i] =  
                pixel.depth = pixel.sqDepth = 0.f;

        }
        pixel.sampleCount = 0;
    }
}

void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->pixelBounds;
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);

        for (int i = 0; i < 3; i++) {
            mergePixel.Lrgb[i] += tilePixel.Lrgb[i];
            mergePixel.sqLrgb[i] += tilePixel.sqLrgb[i];
        }
        mergePixel.sampleCount += tilePixel.sampleCount;
    }
}

void Film::SetImage(const Spectrum *img) const {
}

void Film::AddSplat(const Point2f &p, Spectrum v) {
}

void Film::Update(bool final) {
    for (Point2i p: croppedPixelBounds) {
            auto &pixelInfo = GetPixel(p);
            float invSampleCount = 1.f/(float)pixelInfo.sampleCount;
            float invSampleCount_1 = 1.f/((float)pixelInfo.sampleCount-1.f);
            Color colSum = Color(pixelInfo.Lrgb);
            Color sqColSum = Color(pixelInfo.sqLrgb);
            Color colMean = colSum*invSampleCount;
            Color colVar = (sqColSum - colSum*colMean) * 
                           invSampleCount_1 * invSampleCount;

            Color norSum = Color(pixelInfo.normal);
            Color sqNorSum = Color(pixelInfo.sqNormal);
            Color norMean = norSum*invSampleCount;
            Color norVar = (sqNorSum - norSum*norMean) *
                           invSampleCount_1;
            
            Color rhoSum = Color(pixelInfo.rho);
            Color sqRhoSum = Color(pixelInfo.sqRho);
            Color rhoMean = rhoSum*invSampleCount;
            Color rhoVar = (sqRhoSum - rhoSum*rhoMean) *
                           invSampleCount_1;
            
            float depthSum = pixelInfo.depth;
            float sqDepthSum = pixelInfo.sqDepth;
            float depthMean = depthSum * invSampleCount;
            float depthVar = (sqDepthSum - depthSum*depthMean) *
                             invSampleCount_1;
            
            int x = p.x, y = p.y;
            colImg(x, y) = colMean;
            varImg(x, y) = colVar;
            norImg(x, y) = norMean;
            norVarImg(x, y) = norVar;
            rhoImg(x, y) = rhoMean;            
            rhoVarImg(x, y) = rhoVar;
            depthImg(x, y) = depthMean;
            depthVarImg(x, y) = depthVar;

            Feature feature, featureVar;
            feature[0] = norMean[0];
            feature[1] = norMean[1];
            feature[2] = norMean[2];
            // feature[3] = rhoMean[0];
            // feature[4] = rhoMean[1];
            // feature[5] = rhoMean[2];
            // feature[6] = depthMean;
            featureVar[0] = norVar[0];
            featureVar[1] = norVar[1];
            featureVar[2] = norVar[2];
            // featureVar[3] = rhoVar[0];
            // featureVar[4] = rhoVar[1];
            // featureVar[5] = rhoVar[2];
            // featureVar[6] = depthVar;
            
            featureImg(x, y) = feature;
            featureVarImg(x, y) = featureVar;
    }
    TwoDArray<Color> rColImg = colImg;
    rFilter.Apply(rColImg);
    TwoDArray<Color> rVarImg = varImg;
    rFilter.Apply(rVarImg);


    vector<float> sigma = final ? finalParams : interParams;
    Feature sigmaF;
    sigmaF[0] = sigmaF[1] = sigmaF[2] = sigmaN;

    vector<TwoDArray<Color> > fltArray;
    vector<TwoDArray<float> > mseArray;
    vector<TwoDArray<float> > priArray;
    vector<TwoDArray<float> > wSumArray;
    vector<TwoDArray<float> > fltMseArray;
    vector<TwoDArray<float> > fltPriArray;

    int xPixelCount = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
    int yPixelCount = croppedPixelBounds.pMax.y - croppedPixelBounds.pMin.y;
    for(size_t i = 0; i < sigma.size(); i++) {
        fltArray.push_back(TwoDArray<Color>(xPixelCount, yPixelCount));
        mseArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
        priArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
        wSumArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
        fltMseArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));
        fltPriArray.push_back(TwoDArray<float>(xPixelCount, yPixelCount));

        CrossBilateralFilter cbFilter(sigma[i], 0, sigmaF, xPixelCount, yPixelCount); 
        TwoDArray<Color> flt(xPixelCount, yPixelCount);
        TwoDArray<float> mse(xPixelCount, yPixelCount);            
        TwoDArray<float> pri(xPixelCount, yPixelCount);
        cbFilter.Apply(colImg, featureImg, featureVarImg, rColImg, varImg, rVarImg, flt, mse, pri);
        mseArray[i] = mse;
        priArray[i] = pri;
        fltArray[i] = flt;
    }

    CrossBilateralFilter mseFilter(final ? finalMseSigma : interMseSigma, 0.f, sigmaF, xPixelCount, yPixelCount); 
    mseFilter.Apply(mseArray, priArray, featureImg, featureVarImg, fltMseArray, fltPriArray);

    minMseImg = numeric_limits<float>::infinity();   
    for(size_t i = 0; i < sigma.size(); i++) {
        for(int y = 0; y < yPixelCount; y++)
            for(int x = 0; x < xPixelCount; x++) {
                float error = fltMseArray[i](x, y);                
                float pri = fltPriArray[i](x, y);
                if(error < minMseImg(x, y)) {
                    Color c = fltArray[i](x, y);
                    adaptImg(x, y) = max(pri, 0.f) / (float)(1.f + GetPixel(Point2i(x, y)).sampleCount);
                    minMseImg(x, y) = error;
                    fltImg(x, y) = c;
                    sigmaImg(x, y) = Color((float)i/(float)sigma.size());
            }
        }
    }
}

void Film::GetAdaptPixels(float avgSpp, vector<vector<int> > &pixOff, vector<vector<int> > &pixSmp) {
    Update(false);

    // Clear pixels
    vector<vector<int> >().swap(pixOff);
    vector<vector<int> >().swap(pixSmp);

    int xPixelCount = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
    int yPixelCount = croppedPixelBounds.pMax.y - croppedPixelBounds.pMin.y;

    // Fill offsets
    pixOff.resize(yPixelCount);
    for(int y = 0; y < yPixelCount; y++) {
        pixOff[y].resize(xPixelCount);
        for(int x = 0; x < xPixelCount; x++) 
            pixOff[y][x] = GetPixel(Point2i(x, y)).sampleCount;        
    }
    
    long double totalSamples = (long double)xPixelCount*
                             (long double)yPixelCount*
                             (long double)avgSpp;
    
    long double probSum = 0.0L;
    for(int y = 0; y < yPixelCount; y++)
        for(int x = 0; x < xPixelCount; x++) {
            probSum += (long double)adaptImg(x, y);
        }
    long double invProbSum = 1.0L/probSum;

    pixSmp.resize(yPixelCount);
    for(int y = 0; y < yPixelCount; y++) {
        pixSmp[y].resize(xPixelCount);
        for(int x = 0; x < xPixelCount; x++) {
            pixSmp[y][x] = 
                max((int)ceil(totalSamples * 
                    (long double)adaptImg(x, y) * invProbSum), 1);
        }
    }
}

void Film::WriteImage(Float splatScale) {
    Update(true);
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
        Pixel &pixel = GetPixel(p);
        rgb[3*offset] = fltImg(p.x, p.y)[0];
        rgb[3*offset + 1] = fltImg(p.x, p.y)[1];
        rgb[3*offset + 2] = fltImg(p.x, p.y)[2];
        ++offset;
    }

    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
    pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    std::string filename;
    if (PbrtOptions.imageFile != "") {
        filename = PbrtOptions.imageFile;
        std::string paramsFilename = params.FindOneString("filename", "");
        if (paramsFilename != "")
            Warning(
                "Output filename supplied on command line, \"%s\" is overriding "
                "filename provided in scene description file, \"%s\".",
                PbrtOptions.imageFile.c_str(), paramsFilename.c_str());
    } else
        filename = params.FindOneString("filename", "pbrt.exr");

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop;
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);
    else
        crop = Bounds2f(Point2f(Clamp(PbrtOptions.cropWindow[0][0], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][0], 0, 1)),
                        Point2f(Clamp(PbrtOptions.cropWindow[0][1], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][1], 0, 1)));

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                   Infinity);

    int nInterParams = 0;
    const float *interParams = params.FindFloat("interparams", &nInterParams);
    vector<float> interParamsV;
    if(nInterParams == 0) {
        interParamsV.push_back(0.f);
    } else {
        for(int i = 0; i < nInterParams; i++)
            interParamsV.push_back(interParams[i]);
    }
    int nFinalParams = 0;
    const float *finalParams = params.FindFloat("finalparams", &nFinalParams);
    vector<float> finalParamsV;
    if(nFinalParams == 0) {
        finalParamsV.push_back(0.f);
    } else {
        for(int i = 0; i < nFinalParams; i++)
            finalParamsV.push_back(finalParams[i]);
    }

    float sigmaN = params.FindOneFloat("sigman", 0.8f);
    float sigmaR = params.FindOneFloat("sigmar", 0.25f);
    float sigmaD = params.FindOneFloat("sigmad", 0.6f);
    float interMseSigma = params.FindOneFloat("intermsesigma", 4.f);
    float finalMseSigma = params.FindOneFloat("finalmsesigma", 8.f);

    return new Film(Point2i(xres, yres), crop, std::move(filter), diagonal,
                    filename, scale, interParamsV, finalParamsV, sigmaN, sigmaR, sigmaD, interMseSigma, finalMseSigma, maxSampleLuminance);
}

}  // namespace pbrt
