cmake_minimum_required(VERSION 3.12)

# set the project name
project(PhotoICE VERSION 4.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# pass version number to the header source file
configure_file(photoion_header.h.in photoion_header.h)

# add the executable
add_executable(PhotoICE photoion_comp.cpp)

# Add the binary folder to the include research path
target_include_directories(PhotoICE PUBLIC
                              "${PROJECT_BINARY_DIR}")
