#include "driver/context.h"
#include "gi/algorithm.h"
#include "gi/random.h"
#include "gi/ray.h"
#include "gi/mesh.h"
#include "gi/material.h"
#include "gi/light.h"
#include "gi/timer.h"
#include "gi/color.h"
#include <atomic>

using namespace std;
using namespace glm;

struct NaivePathtracer : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new NaivePathtracer; }

void NaivePathtracer::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    const Camera& cam = context.cam;
    const Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    const size_t w = fbo.width(), h = fbo.height();

    // trace
    for (uint32_t i = 0; i < samples; ++i) {
        vec3 radiance(0);
        // TODO ASSIGNMENT4
        // - implement a (naive) pathtracer using BRDF sampling
        // - add russian roulette


        
        context.fbo.add_sample(x, y, radiance);
    }
}
