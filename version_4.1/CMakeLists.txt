#CMake config
cmake_minimum_required(VERSION 3.12)
set(CMAKE_VERBOSE_MAKEFILE ON)

# set the project name
project(PhotoICE VERSION 4.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler config
set(CMAKE_CXX_COMPILER icpc)
set (CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -v -g -Wall -Wextra -Wpedantic -Werror -Wno-int-in-bool-context")

# require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
       message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
    endif()

# pass version number to the header source file
configure_file(photoion_header.h.in photoion_header.h)

# add the MathFunctions library directory
add_subdirectory(MathFunctions)
add_subdirectory(FilesUtils)
add_subdirectory(Utilities)
add_subdirectory(Tests)


add_library(GSL STATIC IMPORTED)
set_property(TARGET GSL PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgsl.so)
#find_library(GSL_LIBLOC gsl PATHS ${GSL_ROOT} PATH_SUFFIXES lib include)

# add the executables
add_executable(PhotoICE photoion_comp.cpp)


# Add the binary folder to the include research path
target_include_directories(PhotoICE PUBLIC
                              "${PROJECT_BINARY_DIR}"
                              "${PROJECT_SOURCE_DIR}/MathFunctions"
                              "${PROJECT_SOURCE_DIR}/FilesUtils"
                              "${PROJECT_SOURCE_DIR}/Utilities")

target_link_libraries(PhotoICE PUBLIC MathFunctions FilesUtils Utilities)

install(TARGETS PhotoICE DESTINATION bin)
install(FILES "${PROJECT_BINARY_DIR}/Photoion_header.h"
DESTINATION include)
