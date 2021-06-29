#include "driver/context.h"
#include "gi/algorithm.h"
#include "gi/random.h"
#include "gi/ray.h"
#include "gi/mesh.h"
#include "gi/material.h"
#include "gi/light.h"
#include "gi/timer.h"
#include "gi/color.h"
#include "glm/fwd.hpp"
#include <atomic>

using namespace std;
using namespace glm;

struct Pathtracer : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new Pathtracer; }


// radiance ret
glm::vec3 tracePath(Ray& ray, Context& context, int N, float throughput_acc){
    // intersect main ray with scene
    const SurfaceInteraction hit = context.scene.intersect(ray);
    const Scene& scene = context.scene;

    // check if a hit was found
    if (hit.valid) {
        if (hit.is_light() && dot(normalize(hit.N), normalize(-ray.dir)) > 0.f) {// direct light source hit
            return hit.Le();
        }else { // surface hit
            // indirect Light
            auto [brdf, w_i, pdf] = hit.sample(-normalize(ray.dir), RNG::uniform<vec2>());
            vec3 indirect_radiance = vec3(0);
            if (pdf != 0.f){ // there can be no light from here (wrong hemisshere check failed)
                auto newRay = hit.spawn_ray(w_i);
                auto throughput = throughput_acc * (brdf.x + brdf.y + brdf.z) / 3.f * fmaxf(0.f, dot(normalize(hit.N), normalize(newRay.dir))) / pdf;
                auto terminationP = throughput / (1-context.RR_THRESHOLD);
                if(N < context.RR_MIN_PATH_LENGTH || terminationP > context.RR_THRESHOLD){
                    auto Li = tracePath(newRay, context, N+1, throughput);
                    indirect_radiance = brdf * Li * fmaxf(0.f, dot(normalize(hit.N), normalize(newRay.dir))) / pdf;
                    assert(!std::isnan(indirect_radiance.x));
                }else{
                    pdf = 0.f;
                }

            }

            // direct Light
            vec3 direct_radiance = vec3(0);
            const auto [light, pdf_light_source] = scene.sample_light_source(RNG::uniform<float>());
            auto [light_Li, shadow_ray, pdf_light_sample] = light->sample_Li(hit, RNG::uniform<vec2>());
            const float light_pdf = pdf_light_source * pdf_light_sample;
            if (light_pdf > 0.f && !scene.occluded(shadow_ray)){
                direct_radiance = hit.brdf(-ray.dir, shadow_ray.dir) * light_Li * fmaxf(0.f, dot(hit.N, shadow_ray.dir)) / light_pdf;
            }

            auto poor_mans_mip = [](auto pdf_l, auto l, auto pdf_r, auto r){
                if(std::isnan(pdf_l)) return r;
                if(std::isnan(pdf_r)) return l;
                if(pdf_r == 0.f && pdf_l == 0.f) return glm::vec3(0);
                return balance_heuristic(pdf_l, pdf_r) * l + balance_heuristic(pdf_r, pdf_l) * r; 
            };

            //auto radiance = poor_mans_mip(pdf, indirect_radiance, light_pdf, direct_radiance);

            // a bit unsure about this, because the pseudo code is kinda different
            auto radiance = indirect_radiance + direct_radiance; 


            return radiance;
        }
    } else {// ray esacped the scene
        return context.scene.Le(ray);
    }
    return glm::vec3(0);
}

void Pathtracer::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    const Camera& cam = context.cam;
    const Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    const size_t w = fbo.width(), h = fbo.height();
    const auto MAX_PATH_LENGTH = context.MAX_CAM_PATH_LENGTH;


    // trace
    for (uint32_t i = 0; i < samples; ++i) {
        // TODO ASSIGNMENT4
        // - implement a pathtracer using next event estimation
        // - add russian roulette
        // - (optional, bonus) add multiple importance sampling
        Ray ray = cam.view_ray(x, y, w, h, RNG::uniform<vec2>());

        auto radiance = tracePath(ray, context, MAX_PATH_LENGTH, 1.f);

        // add radiance (exitance) to framebuffer
        fbo.add_sample(x, y, radiance);
    }
}
