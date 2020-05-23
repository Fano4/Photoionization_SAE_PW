add_library(gsl SHARED IMPORTED)
set_property(TARGET gsl PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgsl.so)

add_library(gslcblas SHARED IMPORTED)
set_property(TARGET gslcblas PROPERTY IMPORTED_LOCATION $ENV{GSL_ROOT}/lib/libgslcblas.so)

add_library(flint SHARED IMPORTED)
set_property(TARGET flint PROPERTY IMPORTED_LOCATION $ENV{FLINT_ROOT}/lib/libflint.so)

add_library(arb SHARED IMPORTED)
set_property(TARGET arb PROPERTY IMPORTED_LOCATION $ENV{ARB_ROOT}/lib/libarb.so)

add_library(MathFunctions angular_int_aux.cpp Molec_integrals.cpp spherical_harmonics.cpp)
target_include_directories(MathFunctions   SYSTEM    PRIVATE
                                 "$ENV{GSL_ROOT}/include"
                                 "$ENV{MPFR_ROOT}/include"
                                 "$ENV{FLINT_ROOT}/include"
                                 "$ENV{ARB_ROOT}/include")



target_link_libraries(MathFunctions PRIVATE gsl)
target_link_libraries(MathFunctions PRIVATE gslcblas)
target_link_libraries(MathFunctions PRIVATE arb)
target_link_libraries(MathFunctions PRIVATE flint)
install(TARGETS MathFunctions DESTINATION lib)
install(FILES MathFunctions.h DESTINATION include)
