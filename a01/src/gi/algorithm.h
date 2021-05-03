#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <string>

// forward declarations
class Context;
namespace json11 { class Json; }

// ---------------------------------------------------------------------------------------
// XXX implement these in your source file to enable dynamic linking to the module at runtime:
//extern "C" Algorithm* create_algorithm() { return new MyAlgorithm; }
// ---------------------------------------------------------------------------------------

/**
 * @brief Base class for algorithm modules, providing the shared interface
 * @note implement this in your algorithm's source file to enable dynamic linking of the module at runtime:
 * extern "C" Algorithm* create_algorithm() { return new MyAlgorithm; }
 */
class Algorithm {
public:
    virtual ~Algorithm() {}

    /**
     * @brief Read json11::Json config provided by the user, e.g. for algorithm specific parameters
     *
     * @param cfg json11::Json config file to read from
     */
    virtual void read_config(const json11::Json& cfg) {}

    /**
     * @brief Preprocessing step done once before each call to render(), e.g. to build a photon map
     *
     * @param context Context reference, providing the context
     */
    virtual void init(Context& context) {}

    /**
     * @brief Actual render call, responsible for filling the Framebuffer
     *
     * @param context reference to the Context
     */
    virtual void sample_pixel(Context& context, uint32_t x, uint32_t y, uint32_t samples) = 0;
};

// ---------------------------------------------------------------------------------
// module loader helper functions

#if defined(__unix__)
    #include <dlfcn.h>
#elif defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
#else
    #error Unsupported platform!;
#endif

inline void die(const char* msg) {
    fprintf(stderr, "ERROR: %s:%i: %s\n", __FILE__, __LINE__, msg);
    exit(1);
}

inline void* open_module(const char* filepath) {
#if defined(__unix__)
    void* handle = dlopen(filepath, RTLD_NOW);
    if (!handle) die(dlerror());
#else
    void* handle = (void*)LoadLibrary(filepath);
    if (!handle) die((std::string("Unable to load library: ") + filepath).c_str());
#endif
    return handle;
}

inline void close_module(void* handle) {
#if defined(__unix__)
    dlclose(handle);
#else
    FreeLibrary((HINSTANCE)handle);
#endif
}

inline void* resolve_symbol(void* handle, const char* symbol_name) {
#if defined(__unix__)
    void* symbol = dlsym(handle, symbol_name);
    if (!symbol) die(dlerror());
#else
    void* symbol = (void*)GetProcAddress((HINSTANCE)handle, symbol_name);
	if (!symbol) die((std::string("Unable to load symbol: ") + symbol_name).c_str());
#endif
    return symbol;
}
