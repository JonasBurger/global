project(gi)
cmake_minimum_required(VERSION 3.10)

if (NOT UNIX AND NOT WIN32)
	message(FATAL_ERROR "Unsupported platform!")
endif()

# ---------------------------------------------------------------------
# path management

# look for "dep" folder in root or parent directories
set(DEPDIR "${CMAKE_SOURCE_DIR}/dep")
if (NOT EXISTS ${DEPDIR})
    set(DEPDIR "${CMAKE_SOURCE_DIR}/../dep")
elseif (NOT EXISTS ${DEPDIR})
    message(FATAL_ERROR "Dependency directory not found, please download and unzip from StudOn!")
endif()

# look for "data" folder in root or parent directories
set(DATADIR "${CMAKE_SOURCE_DIR}/data")
if (NOT EXISTS ${DATADIR})
    set(DATADIR "${CMAKE_SOURCE_DIR}/../data")
elseif (NOT EXISTS ${DATADIR})
    message(FATAL_ERROR "Data directory not found, please download and unzip from StudOn!")
endif()

# output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

message(STATUS "BINDIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "DEPDIR: ${DEPDIR}")
message(STATUS "DATADIR: ${DATADIR}")
message(STATUS "CONFDIR: ${CMAKE_SOURCE_DIR}/configs")

# ---------------------------------------------------------------------
# dependencies

# setup for included dependencies
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)
set(CMAKE_LINK_DIRECTORIES_BEFORE ON)

include_directories(${DEPDIR}/include)

if (UNIX)
    set(CMAKE_BUILD_RPATH ${DEPDIR}/lin/lib)
elseif(WIN32)
    set(CMAKE_BUILD_RPATH ${DEPDIR}/win/bin)
endif()

# ---------------------------------------------------------------------
# options

option(STAT_COUNTERS "Enable statistic counters." OFF)
option(OPENIMAGEDENOISE "Enable OpenImageDenoise." ON)

if(STAT_COUNTERS)
    add_definitions(-DSTATS)
else()
    add_definitions(-USTATS)
endif()

if(OPENIMAGEDENOISE)
    add_definitions(-DWITH_OIDN)
else()
    add_definitions(-UWITH_OIDN)
endif()

# ---------------------------------------------------------------------
# compiler options

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

add_definitions(-DGI_BIN_DIR="${CMAKE_CURRENT_BINARY_DIR}")
add_definitions(-DGI_DEP_DIR="${DEPDIR}")
add_definitions(-DGI_DATA_DIR="${DATADIR}")
add_definitions(-DGI_CONF_DIR="${CMAKE_SOURCE_DIR}/configs")

if(UNIX)
    add_compile_options(-pthread -fopenmp -march=native -ffast-math -flto)
    add_compile_options(-Wall -Wno-unused-parameter -Wno-unused-variable -Wno-strict-aliasing -Wno-class-memaccess -Wno-misleading-indentation)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON) # somehow breaks auto dll export on windows
else()
    add_definitions(-DNOMINMAX)
    add_definitions(-D_USE_MATH_DEFINES)
    add_definitions(-D_SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING) # wtf
endif()

include(GenerateExportHeader)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# ---------------------------------------------------------------------
# traverse source tree

include_directories(src)
add_subdirectory(src)
