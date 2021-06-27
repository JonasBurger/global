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
    const auto MAX_PATH_LENGHT = context.MAX_CAM_PATH_LENGTH;

    // trace
    for (uint32_t i = 0; i < samples; ++i) {
        vec3 radiance(0);
        // TODO ASSIGNMENT4
        // - implement a (naive) pathtracer using BRDF sampling
        // - add russian roulette

        // setup a view ray
        Ray ray = cam.view_ray(x, y, w, h, RNG::uniform<vec2>(), RNG::uniform<vec2>());
        float fAcc = 1.f;

        for(int i = 0; i < MAX_PATH_LENGHT; ++i){

            // intersect main ray with scene
            const SurfaceInteraction hit = scene.intersect(ray);
            // check if a hit was found
            if (hit.valid) {
                

                if (hit.is_light()) {// direct light source hit
                    radiance += fAcc * hit.Le();
                    break;
                }else { // surface hit
                    auto oldRayDir = ray.dir;
                    const auto [brdf, w_i, pdf] = hit.sample(-ray.dir, RNG::uniform<vec2>());
                    ray.org = hit.P;
                    ray.dir = w_i;
                    //radiance += fAcc * hit.albedo();
                    fAcc *= pdf;
                }
            } else {// ray esacped the scene
                radiance += fAcc * scene.Le(ray);
                break;
            }
        }

        // add radiance (exitance) to framebuffer
        context.fbo.add_sample(x, y, radiance);
    }
}


