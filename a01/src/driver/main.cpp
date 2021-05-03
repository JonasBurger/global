#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#include "context.h"
#include "gi/random.h"
#include "gi/distribution.h"

// ---------------------------------------------------------------------------------
// MAIN

int main(int argc, char** argv) {
    if (!(argc == 2 || argc == 3)) {
        std::cout << "Usage: " << argv[0] << " <path_to_algorithm_module> [path_to_config_file]\n";
        exit(1);
    }

    plot_all_samplers2D();

    // init context and load algorithm/config
    Context ctx;
    ctx.load_algorithm(argv[1]);
    if (argc == 3)
        ctx.load_config(argv[2]);
		
    // run
    ctx.run();
    return 0;
}
