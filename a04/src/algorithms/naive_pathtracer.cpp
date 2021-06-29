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

struct NaivePathtracer : public Algorithm {
    void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples);
};
extern "C" Algorithm* create_algorithm() { return new NaivePathtracer; }

// radiance ret
glm::vec3 tracePath(Ray& ray, Context& context, int N){
    if(N==0){return glm::vec3(0);}
    // intersect main ray with scene
    const SurfaceInteraction hit = context.scene.intersect(ray);
    // check if a hit was found
    if (hit.valid) {
        if (hit.is_light() && dot(normalize(hit.N), normalize(-ray.dir)) > 0.f) {// direct light source hit
            return hit.Le();
        }else { // surface hit
            const auto [brdf, w_i, pdf] = hit.sample(-normalize(ray.dir), RNG::uniform<vec2>());
            if (pdf == 0.f){
                return glm::vec3(0); // there can be no light from here (wrong hemisshere check failed)
            }
            auto newRay = hit.spawn_ray(w_i);
            auto Li = tracePath(newRay, context, N-1);
            auto radiance = brdf * Li * fmaxf(0.f, dot(normalize(hit.N), normalize(newRay.dir))) / pdf; // 
            assert(!std::isnan(radiance.x));
            return radiance;
            //radiance += fAcc * hit.albedo();
        }
    } else {// ray esacped the scene
        return context.scene.Le(ray);
    }
    return glm::vec3(0);
}


void NaivePathtracer::sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) {
    // shortcuts
    const Camera& cam = context.cam;
    const Scene& scene = context.scene;
    Framebuffer& fbo = context.fbo;
    const size_t w = fbo.width(), h = fbo.height();
    const auto MAX_PATH_LENGHT = context.MAX_CAM_PATH_LENGTH;

    // trace
    for (uint32_t i = 0; i < samples; ++i) {
        //vec3 radiance(0);
        // TODO ASSIGNMENT4
        // - implement a (naive) pathtracer using BRDF sampling
        // - add russian roulette

        // setup a view ray
        Ray ray = cam.view_ray(x, y, w, h, RNG::uniform<vec2>());

        auto radiance = tracePath(ray, context, MAX_PATH_LENGHT);

        // add radiance (exitance) to framebuffer
        fbo.add_sample(x, y, radiance);
    }
}



