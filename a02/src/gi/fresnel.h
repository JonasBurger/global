#pragma once

#include "surface.h"

// -------------------------------------------------------------------
// Approximations

inline float fresnel_schlick(float cos_i, float index_of_refraction) {
    // TODO ASSIGNMENT2
    // implement schlick's approximation here
    auto n1 = 1.f; // close to vacuum;
    auto n2 = index_of_refraction;
    auto R0 = pow((n1 - n2)/(n1 + n2), 2);
    auto ROmega = R0 + (1-R0)*pow(1-cos_i, 5);
    return ROmega;
}

// -------------------------------------------------------------------
// Dielectric materials

inline float fresnel_dielectric(float cos_wi, float ior_medium, float ior_material) {
    return fresnel_schlick(cos_wi, ior_material);
}

// -------------------------------------------------------------------
// Conductor materials

inline float fresnel_conductor(float cos_wi, float ior_material, float absorb) {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}
