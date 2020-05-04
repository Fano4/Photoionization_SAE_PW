add_library(gsl SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgsl.so)

add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgslcblas.so)

# add the executables
add_executable(PhotoICE_test photoion_tests.cpp)

target_link_libraries(PhotoICE_test PUBLIC MathFunctions)
#link the directories to the target
target_include_directories(PhotoICE_test PUBLIC
                              "${PROJECT_BINARY_DIR}"
                              "${PROJECT_SOURCE_DIR}/MathFunctions"
                              "${PROJECT_SOURCE_DIR}/Tests"
                              "$ENV{GSL_ROOT}/include")

target_link_libraries(PhotoICE_test PRIVATE gsl)
target_link_libraries(PhotoICE_test PRIVATE gslcblas)
