#include "driver/context.h"
#include "gi/algorithm.h"
#include "gi/sampling.h"
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
        const SurfaceInteraction hit = scene.intersect(ray); // wo
        // check if a hit was found
        if (hit.valid) {
            if (hit.is_light())
                radiance = hit.Le();
            else {
                constexpr bool UNIFORM = true;
                if constexpr (UNIFORM){
                    // TODO ASSIGNMENT2
                    // implement Monte Carlo integration via uniform hemisphere sampling here
                    // - draw a uniform random sample on the hemisphere in tangent space and transform it into world-space
                    // - intersect the ray with the scene and check if you hit a light source
                    // - if a light source was hit, compute the radiance via the given equation
                    //radiance = hit.albedo();
                    const auto sample_tangent_space = uniform_sample_hemisphere(RNG::uniform<vec2>());
                    const auto sample_world_space = hit.to_world(sample_tangent_space);
                    auto sample_ray = Ray(hit.P, sample_world_space); // wi
                    const auto sample_hit = scene.intersect(sample_ray);

                    const float pdf = uniform_hemisphere_pdf();
                    if (pdf > 0.f && sample_hit.light){
                        const auto Le = sample_hit.Le();
                        radiance = hit.brdf(normalize(-ray.dir), normalize(sample_ray.dir)) * Le * fmaxf(0.f, dot(hit.N, sample_ray.dir)) / pdf;
                    }
                    //radiance = hit.albedo();
                }else{
                    // TODO ASSIGNMENT2
                    // implement Monte Carlo integration via BRDF imporance sampling here
                    // - sample the brdf (BRDF::sample) for a outgoing direction instead of uniform sampling of the hemisphere
                    // - intersect the ray with the scene and check if you hit a light source
                    // - if a light source was hit, compute the radiance via the given equation
                    radiance = hit.albedo();
                }
            }
        } else // ray esacped the scene
            radiance = scene.Le(ray);
        // add radiance to framebuffer
        fbo.add_sample(x, y, radiance);
    }
}
