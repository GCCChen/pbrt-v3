
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

#ifndef PBRT_INTEGRATORS_SBFVOLPATH_H
#define PBRT_INTEGRATORS_SBFVOLPATH_H

// integrators/volpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"
#include "sampler.h"
#include "spectrum.h"

namespace pbrt {

// VolPathIntegrator Declarations
class SBFIntegrator : public Integrator {
  public:
    SBFIntegrator(std::shared_ptr<const Camera> &camera,
                  std::shared_ptr<Sampler> sampler,
                  const Bounds2i &pixelBounds, int initSample = 8, float adaptiveSample = 24.f, int maxSample = 1024, int adaptiveIteration =  1)
        :sampler(sampler), camera(camera), pixelBounds(pixelBounds), initSample(initSample), adaptiveSample(adaptiveSample), maxSample(maxSample), adaptiveIteration(adaptiveIteration) {}
    void Render(const Scene &scene);
    virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
    virtual Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth = 0) const = 0;

  private:
    std::shared_ptr<Sampler> sampler;
    std::shared_ptr<const Camera> camera;
    const Bounds2i pixelBounds;
    int initSample, maxSample, adaptiveIteration;
    float adaptiveSample;

};

SBFIntegrator *CreateSBFIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_SBFVOLPATH_H
