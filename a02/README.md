GI - Framework
=======

Global Illumination (GI) framework for use in teaching with focus on GI algorithms instead of BVH construction/traversal.

## Build (Linux)

Install following dependencies via your package manager (all pre-installed in the CIP):

    * Linux C++17 development environment with OpenMP (build-essential, libgomp1, etc)
    * recent GCC with support for std::filesystem
    * OpenGL, GLEW (libglew-dev)
    * GLFW (libglfw3-dev)
    * CMake (cmake)

Download and unpack the `dependencies_linux.zip` from StudOn in the same directory as where you unpacked the assignment `a##.zip` or in the root directory of the source tree.
To build with defaults and out of source, type: `mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -Wno-dev && make` in the root directory.
Use `make -jN` to speed up the build process, where `N` is the amount of threads you want to start.

## Build (Windows)

Install Visual Studio 2019 (https://visualstudio.microsoft.com/de/downloads/) with the C++-packages MSVC and CMake-tools for Windows.
Download and unpack the `dependencies_windows.zip` from StudOn in the same directory as where you unpacked the assignment `a##.zip` or in the root directory of the source tree.
In Visual Studio, open the assignment folder a## (File -> Open -> Folder).
CMake configure and generate should automatically run in the console.
After successful configuration, select the build target `driver` as startup item and hit `Build All`.
The algorithm modules may not be (re-)built automatically when you build the executable, so remember to `Build All`.
Note that the default configuration is Debug, which is several magnitudes slower, so you might want to select Release for rendering.
Additionally, library names may be appended a 'd' in Debug mode, thus you might need to specify `<name>d.algo` instead of `<name>.algo` on the command line in this case.

Alternatively, you can use the program `cmake-gui` to configure the build via selecting the assignment folder a##, a build folder (e.g. a##/build) and clicking on `Configure`, `Generate` and finally `Open Project`.

## Scene Data

Download and unpack the `data.zip` from StudOn in the same directory as where you unpacked the assignment `a##.zip` or in the root directory of the source tree.

## Running

To run, execute the `driver` executable and provide a path to an algorithm module, and optionally JSON configuration, as arguments.
Semantics are: `driver <algorithm module> [configuration file]`.
Note that, if no OpenGL context is available (e.g. when connected to the CIP via SSH), rendering is still possible, albeit without a preview.

## Preview Controls

Use the keys `W A S D R F` to move the camera and drag with the left mouse button pressed to rotate the camera.
You can adjust the camera movement speed via scrolling with the mouse wheel.
Use right mouse click to set the focus plane for depth of field, if the auto focus is disabled.
Press the spacebar to restart rendering and ESC to quit.
See various GUI options to modify the scene and parameters.

## Configs

All options and settings can optionally be imported and exported from/to a single JSON file.
Some example configuration files are located in the `./configs` directory.

## Issues / Suggestions / Feedback

Please mail to <nikolai.hofmann@fau.de>.
