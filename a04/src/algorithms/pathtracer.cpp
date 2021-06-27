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

struct Pathtracer : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new Pathtracer; }

void Pathtracer::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    const Camera& cam = context.cam;
    const Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    const size_t w = fbo.width(), h = fbo.height();

    // trace
    for (uint32_t i = 0; i < samples; ++i) {
        vec3 radiance(0);
        // TODO ASSIGNMENT4
        // - implement a pathtracer using next event estimation
        // - add russian roulette
        // - (optional, bonus) add multiple importance sampling
        context.fbo.add_sample(x, y, radiance);
    }
}
