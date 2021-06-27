#include "brdf.h"
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
    return hit.albedo() * INVPI;
}

std::tuple<vec3, vec3, float> LambertianReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const vec3 w_i = align(hit.N, cosine_sample_hemisphere(sample));
    if (!same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
    return { eval(hit, w_o, w_i), w_i, pdf(hit, w_o, w_i) };
}

float LambertianReflection::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    return abs(dot(hit.N, w_i)) * INVPI;
}

// ----------------------------------------------------------------------------------------------
// Diffuse lambertian transmission

vec3 LambertianTransmission::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return hit.albedo() * INVPI;
}

std::tuple<vec3, vec3, float> LambertianTransmission::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const vec3 w_i = -align(hit.N, cosine_sample_hemisphere(sample));
    if (same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
    return { eval(hit, w_o, w_i), w_i, pdf(hit, w_o, w_i) };
}

float LambertianTransmission::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    return abs(dot(hit.N, w_i)) * INVPI;
}

// ----------------------------------------------------------------------------------------------
// Perfect specular reflection

vec3 SpecularReflection::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return vec3(0);
}

std::tuple<vec3, vec3, float> SpecularReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const vec3 w_i = reflect(-w_o, hit.N);
    if (!same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
    const float cos_i = abs_cos_theta(hit.N, w_i);
    const vec3 brdf = fresnel_dielectric(cos_i, 1.f, hit.mat->ior) * hit.albedo() / cos_i;
    return { brdf, w_i, 1.f };
}

float SpecularReflection::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return 0.f;
}

// ----------------------------------------------------------------------------------------------
// Perfect specular transmission

vec3 SpecularTransmission::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return vec3(0);
}

std::tuple<vec3, vec3, float> SpecularTransmission::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const auto [valid, w_i] = refract(-w_o, hit.N, hit.mat->ior);
    if (!valid || same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
    const float cos_i = abs_cos_theta(hit.N, w_i);
    const vec3 brdf = (1.f - fresnel_dielectric(cos_i, 1.f, hit.mat->ior)) * hit.albedo() / cos_i;
    return { brdf, w_i, 1.f };
}

float SpecularTransmission::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return 0.f;
}

// ----------------------------------------------------------------------------------------------
// Specular fresnel

vec3 SpecularFresnel::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    return vec3(0);
}

std::tuple<vec3, vec3, float> SpecularFresnel::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    const float F = fresnel_dielectric(abs(dot(hit.N, w_o)), 1.f, hit.mat->ior);
    if (sample.x < F) {
        const vec3 w_i = reflect(-w_o, hit.N);
        const vec3 brdf = F * hit.albedo() / abs(dot(hit.N, w_i));
        return { brdf, w_i, F };
    } else {
        const auto [valid, w_i] = refract(-w_o, hit.N, hit.mat->ior);
        if (!valid) return { vec3(0), w_i, 0.f }; // TIR
        const vec3 brdf = (1 - F) * hit.albedo() / abs(dot(hit.N, w_i));
        return { brdf, w_i, 1 - F };
    }
}

float SpecularFresnel::pdf(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    return 0.f;
}

// ----------------------------------------------------------------------------------------------
// Phong

vec3 SpecularPhong::eval(const SurfaceInteraction& hit, const vec3& w_o, const vec3& w_i) const {
    if (!same_hemisphere(hit.N, w_i)) return vec3(0);
    const float exponent = Material::exponent_from_roughness(hit.roughness());
    const float NdotH = fmaxf(0.f, dot(hit.N, normalize(w_o + w_i)));
    const float norm_f = (exponent + 1) * INV2PI;
    const float F = fresnel_dielectric(fmaxf(0.f, abs(dot(hit.N, w_i))), 1.f, hit.mat->ior);
    return F * hit.albedo() * powf(NdotH, exponent) * norm_f;
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
    if (NdotH <= 0) return 0.f;
    const float tan2 = tan2_theta(NdotH);
    if (!std::isfinite(tan2)) return 0.f;
    const float a2 = sqr(roughness);
    assert(M_PI * sqr(sqr(NdotH)) * sqr(a2 + tan2) != 0.f);
    return a2 / (M_PI * sqr(sqr(NdotH)) * sqr(a2 + tan2));
}

inline float GGX_G1(const float NdotV, float roughness) {
    if (NdotV <= 0) return 0.f;
    const float tan2 = tan2_theta(NdotV);
    if (!std::isfinite(tan2)) return 0.f;
    assert(1.f + sqrtf(1.f + sqr(roughness) * tan2) != 0.f);
    return 2.f / (1.f + sqrtf(1.f + sqr(roughness) * tan2));
}

vec3 GGX_sample(const vec2& sample, float roughness) {
    assert(sample.x != 0.f);
    const float theta = atanf((roughness * sqrtf(sample.x)) / sqrtf(1.f - sample.x));
    if (!std::isfinite(theta)) return vec3(0, 0, 1);
    const float phi = 2.f * M_PI * sample.y;
    const float sin_t = sinf(theta);
    return vec3(sin_t * cosf(phi), sin_t * sinf(phi), cosf(theta));
}

inline float GGX_pdf(float D, float NdotH, float HdotV) {
    assert(4.f * abs(HdotV) != 0.f);
    return (D * abs(NdotH)) / (4.f * abs(HdotV));
}

// ----------------------------------------------------------------------------------------------
// Microfacet reflection

vec3 MicrofacetReflection::eval(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    const float NdotV = dot(hit.N, w_o);
    const float NdotL = dot(hit.N, w_i);
    if (NdotV <= 0.f || NdotL <= 0.f) return vec3(0);
    const vec3 H = normalize(w_o + w_i);
    const float NdotH = dot(hit.N, H);
    const float HdotL = dot(H, w_i);
    const float roughness = hit.roughness();
    const float F = fresnel_dielectric(HdotL, 1.f, hit.mat->ior);
    const float D = GGX_D(NdotH, roughness);
    const float G = GGX_G1(NdotV, roughness) * GGX_G1(NdotL, roughness);
    assert(abs(NdotV) != 0.f && abs(NdotL) != 0.f);
    const float microfacet = (F * D * G) / (4 * abs(NdotV) * abs(NdotL));
    return coated ? vec3(microfacet) : hit.albedo() * microfacet;
}

std::tuple<vec3, vec3, float> MicrofacetReflection::sample(const SurfaceInteraction& hit, const vec3& w_o, const vec2& sample) const {
    // reflect around sampled macro normal w_h
    const vec3 w_h = align(hit.N, GGX_sample(sample, hit.roughness()));
    const vec3 w_i = reflect(-w_o, w_h);
    if (!same_hemisphere(hit.Ng, w_i)) return { vec3(0), w_i, 0.f };
    const float sample_pdf = pdf(hit, w_o, w_i);
    assert(std::isfinite(sample_pdf));
    return { eval(hit, w_o, w_i), w_i, sample_pdf };
}

float MicrofacetReflection::pdf(const SurfaceInteraction &hit, const vec3 &w_o, const vec3 &w_i) const {
    const vec3 H = normalize(w_o + w_i);
    const float NdotH = dot(hit.N, H);
    const float HdotV = dot(H, w_o);
    const float D = GGX_D(NdotH, hit.roughness());
    return GGX_pdf(D, NdotH, HdotV);
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
    const float NdotV = dot(hit.N, w_o);
    const float NdotL = dot(hit.N, w_i);
    if (NdotV <= 0.f || NdotL <= 0.f) return vec3(0);
    const vec3 H = normalize(w_o + w_i);
    const float NdotH = dot(hit.N, H);
    const float HdotL = dot(H, w_i);
    const float roughness = hit.roughness();
    const float F = fresnel_conductor(HdotL, hit.mat->ior, hit.mat->absorb);
    const float D = GGX_D(NdotH, roughness);
    const float G = GGX_G1(NdotV, roughness) * GGX_G1(NdotL, roughness);
    const float microfacet = (F * D * G) / (4 * abs(NdotV) * abs(NdotL));
    return hit.albedo() * microfacet;
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
