# target name to generate here
set(TARGET driver)

# list source files (via suboptimal lazy glob)
file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

# define target
add_executable(${TARGET} ${SOURCES} ${HEADERS})
target_link_libraries(${TARGET} gi)

# ----------------------------------------------------------
# dependencies

# threads / OpenMP
find_package(Threads REQUIRED)
target_link_libraries(${TARGET} ${CMAKE_THREAD_LIBS_INIT})
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
include_directories(${TARGET} ${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(${TARGET} ${OpenMP_CXX_LIBRARIES})

# OpenGL
set(OpenGL_GL_PREFERENCE "GLVND")
find_package(OpenGL REQUIRED)
include_directories(${TARGET} ${OPENGL_INCLUDE_DIR})
target_link_libraries(${TARGET} ${OPENGL_LIBRARIES})

if (UNIX)
	# GLEW
	find_package(GLEW REQUIRED)
	include_directories(${TARGET} ${GLEW_INCLUDE_DIR})
	target_link_libraries(${TARGET} ${GLEW_LIBRARIES})
	# GLFW3
	find_package(glfw3 REQUIRED)
	target_link_libraries(${TARGET} glfw)
    # std::filesystem
	target_link_libraries(${TARGET} stdc++fs)
elseif (WIN32)
    # GLEW and GLFW
    find_library(GLEW glew32 PATHS ${DEPDIR}/win/lib)
    find_library(GLFW glfw3 PATHS ${DEPDIR}/win/lib)
    message(STATUS "GLEW: ${GLEW}")
    message(STATUS "GLFW: ${GLFW}")
    target_link_libraries(${TARGET} ${GLEW} ${GLFW})
endif()

# dl libraries
target_link_libraries(${TARGET} ${CMAKE_DL_LIBS})

# windows: copy dlls to executable dir
if (WIN32)
	add_custom_command(TARGET ${TARGET} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory "${DEPDIR}/win/bin" $<TARGET_FILE_DIR:${TARGET}>)
endif()
