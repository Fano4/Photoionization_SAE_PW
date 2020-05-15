add_library(gsl SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgsl.so)

add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgslcblas.so)

add_library(MathFunctions angular_int_aux.cpp Molec_integrals.cpp spherical_harmonics.cpp)
target_include_directories(MathFunctions       PRIVATE
                                 "$ENV{GSL_ROOT}/include")


target_link_libraries(MathFunctions PRIVATE gsl)
target_link_libraries(MathFunctions PRIVATE gslcblas)
install(TARGETS MathFunctions DESTINATION lib)
install(FILES MathFunctions.h DESTINATION include)
