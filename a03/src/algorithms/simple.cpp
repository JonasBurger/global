#include "driver/context.h"
#include "gi/algorithm.h"
#include "gi/surface.h"
#include "gi/scene.h"
#include "gi/random.h"
#include "gi/light.h"
#include "gi/ray.h"

using namespace std;
using namespace glm;

struct SimpleRenderer : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new SimpleRenderer; }

void SimpleRenderer::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    const Camera& cam = context.cam;
    const Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    const size_t w = fbo.width(), h = fbo.height();

    for (uint32_t i = 0; i < samples; ++i) {
        vec3 radiance(0);
        // setup a view ray
        Ray ray = cam.view_ray(x, y, w, h, RNG::uniform<vec2>(), RNG::uniform<vec2>());
        // intersect main ray with scene
        const SurfaceInteraction hit = scene.intersect(ray);
        // check if a hit was found
        if (hit.valid) {
            if (hit.is_light()) // direct light source hit
                radiance = hit.Le();
            else { // surface hit -> shading
                const auto [light, pdf_light_source] = scene.sample_light_source(RNG::uniform<float>());
                auto [Li, shadow_ray, pdf_light_sample] = light->sample_Li(hit, RNG::uniform<vec2>());
                const float pdf = pdf_light_source * pdf_light_sample;
                if (pdf > 0.f && !scene.occluded(shadow_ray))
                    radiance = hit.brdf(-ray.dir, shadow_ray.dir) * Li * fmaxf(0.f, dot(hit.N, shadow_ray.dir)) / pdf;
            }
        } else // ray esacped the scene
            radiance = scene.Le(ray);
        // add radiance (exitance) to framebuffer
        fbo.add_sample(x, y, radiance);
    }
}
