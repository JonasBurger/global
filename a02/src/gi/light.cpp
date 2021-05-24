#include "light.h"
#include "ray.h"
#include "scene.h"
#include "texture.h"
#include "buffer.h"
#include "mesh.h"
#include "material.h"
#include "sampling.h"
#include "distribution.h"
#include "timer.h"
#include "color.h"
#include <iostream>

// ------------------------------------------------
// Mesh area light

AreaLight::AreaLight(const Mesh& mesh) : mesh(mesh) {}

std::tuple<glm::vec3, Ray, float> AreaLight::sample_Li(const SurfaceInteraction& hit, const glm::vec2& sample) const {
    assert(sample.x >= 0 && sample.x < 1); assert(sample.y >= 0 && sample.y < 1);
    STAT("sampleLi");
    // sample area light source (triangle mesh)
    const auto [light, sample_pdf] = mesh.sample(sample);
    if (sample_pdf <= 0.f) return { glm::vec3(0.f), Ray(), 0.f };
    glm::vec3 l = light.P - hit.P;
    const float r = length(l);
    l = normalize(l);
    const float cos_t_light = dot(light.N, -l);
    if (cos_t_light <= 0.f) return { glm::vec3(0.f), Ray(), 0.f };
    const float pdf = sample_pdf * (r * r) / (cos_t_light * light.area);
    assert(std::isfinite(pdf));
    return { light.Le(), hit.spawn_ray(l, r), pdf };
}

float AreaLight::pdf_Li(const SurfaceInteraction& light, const Ray& ray) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<glm::vec3, Ray, glm::vec3, float, float> AreaLight::sample_Le(const glm::vec2& sample_pos, const glm::vec2& sample_dir) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<float, float> AreaLight::pdf_Le(const SurfaceInteraction& light, const glm::vec3& dir) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

glm::vec3 AreaLight::power() const {
    return glm::vec3(mesh.mat->emissive_strength * mesh.surface_area() * PI);
}

// ------------------------------------------------
// Sky light

SkyLight::SkyLight() {}

SkyLight::SkyLight(const fs::path& path, const Scene& scene, float intensity) {
    load(path, scene.center, scene.radius, intensity);
}

void SkyLight::load(const fs::path& path, const glm::vec3& scene_center, float scene_radius, float intensity) {
    // init variables
    const fs::path resolved_path = fs::exists(path) ? path : fs::path(GI_DATA_DIR) / path;
    std::cout << "loading: " << path << " (" << resolved_path << ")..." << std::endl;
    tex = std::make_shared<Texture>(resolved_path);
    this->intensity = intensity;
    this->scene_center = scene_center;
    this->scene_radius = scene_radius;
    commit();
}

void SkyLight::commit() {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<glm::vec3, Ray, float> SkyLight::sample_Li(const SurfaceInteraction& hit, const glm::vec2& sample) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

float SkyLight::pdf_Li(const SurfaceInteraction& light, const Ray& ray) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

glm::vec3 SkyLight::Le(const Ray& ray) const {
    assert(tex && distribution);
    return tex->env(ray.dir) * intensity;
}

std::tuple<glm::vec3, Ray, glm::vec3, float, float> SkyLight::sample_Le(const glm::vec2& sample_pos, const glm::vec2& sample_dir) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

std::tuple<float, float> SkyLight::pdf_Le(const SurfaceInteraction& light, const glm::vec3& dir) const {
    throw std::runtime_error("Function not implemented: " + std::string(__FILE__) + ", line: " + std::to_string(__LINE__));
}

glm::vec3 SkyLight::power() const {
    assert(tex && distribution);
    return glm::vec3(PI * scene_radius * scene_radius * intensity * distribution->unit_integral());
}

inline std::string fix_data_path(std::string path) {
    // remove GI_DATA_DIR from path to avoid absolute paths
    if (path.find(GI_DATA_DIR) != std::string::npos)
        path = path.substr(path.find(GI_DATA_DIR) + std::strlen(GI_DATA_DIR) + 1);
    return path;
}

json11::Json SkyLight::to_json() const {
    return json11::Json::object{
        { "envmap", fix_data_path(tex->path().string()) },
        { "intensity", intensity },
        { "scene_center", json11::Json::array{scene_center.x, scene_center.y, scene_center.z} },
        { "scene_radius", scene_radius },
    };
}

void SkyLight::from_json(const json11::Json& cfg) {
    if (cfg.is_object()) {
        json_set_float(cfg, "intensity", intensity);
        json_set_vec3(cfg, "scene_center", scene_center);
        json_set_float(cfg, "scene_radius", scene_radius);
        if (cfg["envmap"].is_string())
            load(cfg["envmap"].string_value(), scene_center, scene_radius, intensity);
        else
            tex = std::make_shared<Texture>(glm::vec3(1));
        commit();
    }
}
