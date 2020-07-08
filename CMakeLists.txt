cmake_minimum_required(VERSION 3.8)
project(delaunator VERSION 0.4.0)

set(ROOT_DIR ${PROJECT_SOURCE_DIR})

#
# Must be first.
#
include(${ROOT_DIR}/cmake/common.cmake NO_POLICY_SCOPE)

include(CTest)
include(FeatureSummary)

include(${D_CMAKE_DIR}/options.cmake)
include(${D_CMAKE_DIR}/gtest.cmake)

configure_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/delaunator_config.hpp.in
    ${CMAKE_CURRENT_LIST_DIR}/include/delaunator_config.hpp
    @ONLY
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${D_BIN_DIR})

if (DELAUNATOR_HEADER_ONLY)
    add_library(delaunator INTERFACE)
    target_include_directories(delaunator INTERFACE include/)
else()
    add_library(delaunator include/delaunator.cpp)
    target_include_directories(delaunator PUBLIC include/)
    d_target_compile_settings(delaunator)
endif()

add_subdirectory(test)

#examples
#add_executable(triangulate-geojson examples/triangulate_geojson.cpp include/delaunator.cpp)
#add_executable(basic examples/basic.cpp include/delaunator.cpp)

#set(CPACK_PROJECT_NAME ${PROJECT_NAME})
#set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
