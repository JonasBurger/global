#include "driver/context.h"
#include "gi/algorithm.h"
#include "gi/surface.h"
#include "gi/scene.h"
#include "gi/random.h"
#include "gi/light.h"
#include "gi/ray.h"

using namespace std;
using namespace glm;

struct BRDFImportance : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new BRDFImportance; }

void BRDFImportance::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    Camera& cam = context.cam;
    Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    size_t w = fbo.width(), h = fbo.height();

    for (uint32_t i = 0; i < samples; ++i) {
        vec3 radiance(0);
        // setup view ray
        Ray ray = cam.view_ray(x, y, w, h, RNG::uniform<vec2>(), RNG::uniform<vec2>());
        // intersect main ray with scene
        const SurfaceInteraction hit = scene.intersect(ray);
        // check if a hit was found
        if (hit.valid) {
            if (hit.is_light())
                radiance = hit.Le();
            else {
                #define UNIFORM
                #if defined(UNIFORM)
                    // TODO ASSIGNMENT2
                    // implement Monte Carlo integration via uniform hemisphere sampling here
                    // - draw a uniform random sample on the hemisphere in tangent space and transform it into world-space
                    // - intersect the ray with the scene and check if you hit a light source
                    // - if a light source was hit, compute the radiance via the given equation
                    radiance = hit.albedo();
                #else
                    // TODO ASSIGNMENT2
                    // implement Monte Carlo integration via BRDF imporance sampling here
                    // - sample the brdf (BRDF::sample) for a outgoing direction instead of uniform sampling of the hemisphere
                    // - intersect the ray with the scene and check if you hit a light source
                    // - if a light source was hit, compute the radiance via the given equation
                    radiance = hit.albedo();
                #endif
            }
        } else // ray esacped the scene
            radiance = scene.Le(ray);
        // add radiance to framebuffer
        fbo.add_sample(x, y, radiance);
    }
}
