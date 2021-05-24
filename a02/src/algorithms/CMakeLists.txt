# list source files (via suboptimal lazy glob)
file(GLOB SOURCES "*.cpp")

# ----------------------------------------------------------
# library targets

foreach(SOURCE ${SOURCES})
    get_filename_component(ALGO ${SOURCE} NAME_WE)
    add_library(${ALGO} SHARED ${SOURCE})
    generate_export_header(${ALGO})
    target_link_libraries(${ALGO} gi)
    if (UNIX)
        target_link_libraries(${ALGO} ${CMAKE_THREAD_LIBS_INIT} ${OpenMP_CXX_LIBRARIES})
    endif()
    set_target_properties(${ALGO} PROPERTIES PREFIX "")
    set_target_properties(${ALGO} PROPERTIES SUFFIX ".algo")
endforeach()
