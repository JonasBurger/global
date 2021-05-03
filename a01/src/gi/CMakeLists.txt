# target name to generate here
set(TARGET gi)

# list source files (via suboptimal lazy glob)
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

# define target
add_library(${TARGET} SHARED ${SOURCES} ${HEADERS})
generate_export_header(${TARGET})

# ----------------------------------------------------------
# dependencies

# threads / OpenMP
find_package(Threads REQUIRED)
target_link_libraries(${TARGET} ${CMAKE_THREAD_LIBS_INIT})
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
include_directories(${TARGET} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(${TARGET} ${OpenMP_CXX_LIBRARIES})

# embree
find_library(EMBREE embree3 PATHS ${DEPDIR}/lin/lib ${DEPDIR}/win/lib NO_DEFAULT_PATH)
message(STATUS "EMBREE: ${EMBREE}")
target_link_libraries(${TARGET} ${EMBREE})

# ASS-imp TODO windows
find_library(ASSIMP assimp PATHS ${DEPDIR}/lin/lib ${DEPDIR}/win/lib NO_DEFAULT_PATH)
message(STATUS "ASSIMP: ${ASSIMP}")
target_link_libraries(${TARGET} ${ASSIMP})

# OIDN
if(OPENIMAGEDENOISE)
    find_library(OIDN OpenImageDenoise PATHS ${DEPDIR}/lin/lib ${DEPDIR}/win/lib NO_DEFAULT_PATH)
    message(STATUS "OIDN: ${OIDN}")
    target_link_libraries(${TARGET} ${OIDN})
endif()
