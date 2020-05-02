#ifndef FRONUC_DYN_M_H
#define FRONUC_DYN_M_H
#define MKL_Complex16 std::complex<double>


//Include standards libraries
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <ctime>
#include "mkl.h"
#include "omp.h"

//include program headers files
#include "prime.hpp"
#include "files_reader.hpp"
//#include "dyson_cube_writer.hpp"
#include "algebra.hpp"
#include "Computation.hpp"
#include "hf5_photoion.hpp"
#include "wignerSymbols-master/src/wignerSymbols-cpp.cpp"
#include "files_reader.cpp"
//#include "dyson_cube_writer.cpp"
#include "algebra.cpp"
#include "angular_int_aux.cpp"
#include "Computation.cpp"
#include "hf5_photoion.cpp"

int spherical_harmonics_translator(std::string basis_func_type,bool component);

#endif
