cmake_minimum_required(VERSION 2.6)
project(complementary_dynamics)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED)

# Add your project files
include_directories("include/")
file(GLOB SRCFILES src/*.cpp)

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::opengl igl::opengl_glfw )

add_executable(complementary_dynamics "main.cpp")
target_link_libraries(complementary_dynamics core igl::core igl::opengl igl::opengl_glfw)
