#include "brdf.h"
#include "glm/geometric.hpp"
#include "surface.h"
#include "fresnel.h"
#include "material.h"
#include "sampling.h"
#include "color.h"
#include <cmath>

using namespace glm;

// ----------------------------------------------------------------------------------------------
// Diffuse lambertian reflection

vec3 LambertianReflection::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    // TODO ASSIGNMENT2
    // evaluate the (normalized!) lambertian diffuse BRDF
    //return vec3(0);
    vec3 CIl = hit.Le();
    vec3 N = hit.N;
    vec3 L = w_o;
    auto ID = dot(L, N)*CIl;
    // todo: normalize?
    return ID;
}

std::tuple<vec3, vec3, float> LambertianReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    // TODO ASSIGNMENT2
    // importance sample and evaluate the lambertian diffuse BRDF
    // set w_i to the sampled (world-space!) direction, pdf to the respective PDF and brdf to the evaluated BRDF
    const vec3 w_i = vec3(0);
    const vec3 brdf = vec3(0);
    const float pdf = 0.f;
    return { brdf, w_i, pdf };
}

float LambertianReflection::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Diffuse lambertian transmission

vec3 LambertianTransmission::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> LambertianTransmission::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float LambertianTransmission::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Perfect specular reflection

vec3 SpecularReflection::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> SpecularReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float SpecularReflection::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Perfect specular transmission

vec3 SpecularTransmission::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> SpecularTransmission::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float SpecularTransmission::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Specular fresnel

vec3 SpecularFresnel::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> SpecularFresnel::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float SpecularFresnel::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Phong

vec3 SpecularPhong::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    // TODO ASSIGNMENT2
    // evaluate the (normalized!) phong BRDF for the given in- and outgoing (world-space) directions
    const float exponent = Material::exponent_from_roughness(hit.roughness());
    const float index_of_refraction = hit.mat->ior;
    return vec3(0);
}

std::tuple<vec3, vec3, float> SpecularPhong::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float SpecularPhong::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Microfacet distribution helper functions

inline float GGX_D(const float NdotH, float roughness) {
    // TODO ASSIGNMENT2 (optional)
    // compute the GGX D term here
    return 0.f;
}

inline float GGX_G1(const float NdotV, float roughness) {
    // TODO ASSIGNMENT2 (optional)
    // compute the GGX G1 term here
    return 0.f;
}

vec3 GGX_sample(const vec2& sample, float roughness) {
    // TODO ASSIGNMENT2 (optional)
    // implement sampling the GGX distribution here
    // note: return a direction in tangent space
    return vec3(0);
}

inline float GGX_pdf(float D, float NdotH, float HdotV) {
    // TODO ASSIGNMENT2 (optional)
    // compute the microfacet PDF here
    return 0.f;
}

// ----------------------------------------------------------------------------------------------
// Microfacet reflection

vec3 MicrofacetReflection::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    // TODO ASSIGNMENT2
    // evaluate the full microfacet BRDF here, optionally relying on the above functions for the D and G1 terms
    // note: use schlick's approximation for the F term
    const float alpha = hit.roughness();
    const float microfacet = 0.f;
    return coated ? vec3(microfacet) : hit.albedo() * microfacet;
}

std::tuple<vec3, vec3, float> MicrofacetReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    // TODO ASSIGNMENT2
    // importance sample and evaluate this microfacet BRDF
    // set w_i to the sampled (world-space!) direction, pdf to the respective PDF and brdf to the evaluated BRDF
    const vec3 w_i = vec3(0);
    const vec3 brdf = vec3(0);
    const float pdf = 0.f;
    return { brdf, w_i, pdf };
}

float MicrofacetReflection::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ------------------------------------------------
// Microfacet transmission

vec3 MicrofacetTransmission::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> MicrofacetTransmission::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float MicrofacetTransmission::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// -------------------------------------------------------------------------------------------
// Layered

vec3 LayeredSurface::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    const float F = fresnel_dielectric(abs(dot(hit.N, w_o)), 1.f, hit.mat->ior);
    return mix(diff.eval(hit, w_o, w_i), spec.eval(hit, w_o, w_i), F);
}

std::tuple<vec3, vec3, float> LayeredSurface::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const float F = fresnel_dielectric(abs(dot(hit.N, w_o)), 1.f, hit.mat->ior);
    vec3 brdf;
    if (sample.x < F) {
        // sample specular
        const vec2 sample_mapped = vec2((F - sample.x) / F, sample.y);
        const auto [specular, w_i, sample_pdf] = spec.sample(hit, w_o, sample_mapped);
        if (!same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
        assert(std::isfinite(sample_pdf));
        return { mix(diff.eval(hit, w_o, w_i), specular, F), w_i, mix(diff.pdf(hit, w_o, w_i), sample_pdf, F) };
    } else {
        // sample diffuse
        const vec2 sample_mapped = vec2((sample.x - F) / (1 - F), sample.y);
        const auto [diffuse, w_i, sample_pdf] = diff.sample(hit, w_o, sample_mapped);
        if (!same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
        assert(std::isfinite(sample_pdf));
        return { mix(diffuse, spec.eval(hit, w_o, w_i), F), w_i, mix(sample_pdf, spec.pdf(hit, w_o, w_i), F) };
    }
}

float LayeredSurface::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    const float F = fresnel_dielectric(abs(dot(hit.N, w_o)), 1.f, hit.mat->ior);
    return mix(diff.pdf(hit, w_o, w_i), spec.pdf(hit, w_o, w_i), F);
}

// ----------------------------------------------------------------------------------------------
// Metal

vec3 MetallicSurface::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

// ----------------------------------------------------------------------------------------------
// Glass

vec3 GlassSurface::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<vec3, vec3, float> GlassSurface::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float GlassSurface::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}
