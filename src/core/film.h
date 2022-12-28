
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_FILM_H
#define PBRT_CORE_FILM_H

// core/film.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "stats.h"
#include "parallel.h"

namespace pbrt {

// FilmTilePixel Declarations
struct FilmTilePixel {
    FilmTilePixel() {
        for(int i = 0; i < 3; i++) {
            Lrgb[i] = sqLrgb[i] =
                normal[i] = sqNormal[i] = 
                rho[i] = sqRho[i] =  
                depth = sqDepth = 0.f;
        }
        sampleCount = 0;
    }
    float Lrgb[3];
    float sqLrgb[3];
    float normal[3];
    float sqNormal[3];
    float rho[3];
    float sqRho[3];
    float depth;
    float sqDepth;
    float weightSum;
    int sampleCount;
};

// Film Declarations
class Film {
  public:
    // Film Public Methods
    Film(const Point2i &resolution, const Bounds2f &cropWindow,
         std::unique_ptr<Filter> filter, Float diagonal,
         const std::string &filename, Float scale,
         Float maxSampleLuminance = Infinity);
    Bounds2i GetSampleBounds() const;
    Bounds2f GetPhysicalExtent() const;
    std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);
    void MergeFilmTile(std::unique_ptr<FilmTile> tile);
    void SetImage(const Spectrum *img) const;
    void AddSplat(const Point2f &p, Spectrum v);
    void WriteImage(Float splatScale = 1);
    void Clear();

    // Film Public Data
    const Point2i fullResolution;
    const Float diagonal;
    std::unique_ptr<Filter> filter;
    const std::string filename;
    Bounds2i croppedPixelBounds;

  private:
    // Film Private Data
    struct Pixel {
        Pixel() {
            for(int i = 0; i < 3; i++) {
                Lrgb[i] = sqLrgb[i] =
                    normal[i] = sqNormal[i] = 
                    rho[i] = sqRho[i] =  
                    depth = sqDepth = 0.f;
            }
            sampleCount = 0;
        }
        float Lrgb[3];
        float sqLrgb[3];
        float normal[3];
        float sqNormal[3];
        float rho[3];
        float sqRho[3];
        float depth;
        float sqDepth;
        float weightSum;
        int sampleCount;
    };
    std::unique_ptr<Pixel[]> pixels;
    static PBRT_CONSTEXPR int filterTableWidth = 16;
    Float filterTable[filterTableWidth * filterTableWidth];
    std::mutex mutex;
    const Float scale;
    const Float maxSampleLuminance;

    // Film Private Methods
    Pixel &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        return pixels[offset];
    }
};

class FilmTile {
  public:
    // FilmTile Public Methods
    FilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
             const Float *filterTable, int filterTableSize,
             Float maxSampleLuminance)
        : pixelBounds(pixelBounds),
          filterRadius(filterRadius),
          invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
          filterTable(filterTable),
          filterTableSize(filterTableSize),
          maxSampleLuminance(maxSampleLuminance) {
        pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
    }
    void AddSample(const Point2f &pFilm, Spectrum L, Float sampleWeight = 1.) {
        ProfilePhase _(Prof::AddFilmSample);
        int x = std::floor(pFilm.x);
        int y = std::floor(pFilm.y);

        FilmTilePixel &pixel = GetPixel(Point2i(x, y));

        float rgb[3];
        L.ToRGB(rgb);
        float rhoRGB[3];

        for(int i = 0; i < 3; i++) {
            pixel.Lrgb[i] += rgb[i];
            pixel.sqLrgb[i] += rgb[i]*rgb[i];
            // lack of information about isec
        }    
        pixel.sampleCount++;
    }
    FilmTilePixel &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }
    const FilmTilePixel &GetPixel(const Point2i &p) const {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }
    Bounds2i GetPixelBounds() const { return pixelBounds; }

  private:
    // FilmTile Private Data
    const Bounds2i pixelBounds;
    const Vector2f filterRadius, invFilterRadius;
    const Float *filterTable;
    const int filterTableSize;
    std::vector<FilmTilePixel> pixels;
    const Float maxSampleLuminance;
    friend class Film;
};

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

}  // namespace pbrt

#endif  // PBRT_CORE_FILM_H
