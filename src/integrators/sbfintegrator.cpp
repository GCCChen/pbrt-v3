
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

// integrators/volpath.cpp*
#include "integrators/sbfintegrator.h"
#include <iostream>
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"
#include "geometry.h"
#include "progressreporter.h"

namespace pbrt {
STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);


void SBFIntegrator::Render(const Scene &scene) {
    Preprocess(scene, *sampler);

    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);
    ProgressReporter reporter(nTiles.x * nTiles.y * adaptiveIteration, "Rendering");

    std::cout << "adapt" <<  adaptiveIteration << std::endl;
    

    {
        for (int iterationId = 0; iterationId < adaptiveIteration; iterationId++) {
            vector<vector<int> > pixOff;
            vector<vector<int> > pixSmp;
            camera->film->GetAdaptPixels(adaptiveSample, pixOff, pixSmp);
            auto renderTask = [&](Point2i tile) {
                // Render section of image corresponding to _tile_

                // Allocate _MemoryArena_ for tile
                MemoryArena arena;

                // Get sampler instance for tile
                int seed = tile.y * nTiles.x + tile.x + (nTiles.x * nTiles.y) * iterationId;
                std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

                // Compute sample bounds for tile
                int x0 = sampleBounds.pMin.x + tile.x * tileSize;
                int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
                int y0 = sampleBounds.pMin.y + tile.y * tileSize;
                int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
                Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
                LOG(INFO) << "Starting image tile " << tileBounds;

                // Get _FilmTile_ for tile
                std::unique_ptr<FilmTile> filmTile =
                    camera->film->GetFilmTile(tileBounds);

                // Loop over pixels in tile to render them
                for (Point2i pixel : tileBounds) {
                    {
                        ProfilePhase pp(Prof::StartPixel);
                        tileSampler->StartPixel(pixel);
                        tileSampler->SetSampleNumber(pixOff[pixel.y][pixel.x]);
                    }

                    if (!InsideExclusive(pixel, pixelBounds))
                        continue;

                    int sampleCount = min(iterationId == 0 ? initSample : pixSmp[pixel.y][pixel.x], maxSample);
                    for (int sampleId = 0; sampleId < sampleCount; sampleId++) {
                        // Initialize _CameraSample_ for current sample
                        CameraSample cameraSample =
                            tileSampler->GetCameraSample(pixel);

                        // Generate camera ray for current sample
                        RayDifferential ray;
                        Float rayWeight =
                            camera->GenerateRayDifferential(cameraSample, &ray);
                        ray.ScaleDifferentials(
                            1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                        ++nCameraRays;

                        // Evaluate radiance along camera ray
                        Spectrum L(0.f);
                        if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);

                        // Issue warning if unexpected radiance value returned
                        if (L.HasNaNs()) {
                            LOG(ERROR) << StringPrintf(
                                "Not-a-number radiance value returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        } else if (L.y() < -1e-5) {
                            LOG(ERROR) << StringPrintf(
                                "Negative luminance value, %f, returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                L.y(), pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        } else if (std::isinf(L.y())) {
                                LOG(ERROR) << StringPrintf(
                                "Infinite luminance value returned "
                                "for pixel (%d, %d), sample %d. Setting to black.",
                                pixel.x, pixel.y,
                                (int)tileSampler->CurrentSampleNumber());
                            L = Spectrum(0.f);
                        }
                        VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                            ray << " -> L = " << L;

                        // Add camera ray's contribution to image
                        SurfaceInteraction isect;
                        scene.Intersect(ray, &isect);
                        isect.shadingN = Faceforward(isect.shading.n, ray.d);
                        filmTile->AddSample(cameraSample.pFilm, L, rayWeight, &isect);

                        // Free _MemoryArena_ memory from computing image sample
                        // value
                        arena.Reset();
                        tileSampler->StartNextSample();
                    }
                }
                LOG(INFO) << "Finished image tile " << tileBounds;

                // Merge image tile into _Film_
                camera->film->MergeFilmTile(std::move(filmTile));
                reporter.Update();
            };
            ParallelFor2D(renderTask, nTiles);
        }
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
    camera->film->WriteImage();
  
}

}  // namespace pbrt
