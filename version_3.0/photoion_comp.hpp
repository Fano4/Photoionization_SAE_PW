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
#include "mkl.h"
#include "omp.h"

//include program headers files
#include "files_reader.hpp"
#include "dyson_cube_writer.hpp"
#include "algebra.hpp"
#include "Computation.hpp"

#include "files_reader.cpp"
#include "dyson_cube_writer.cpp"
#include "algebra.cpp"
#include "Computation.cpp"


#endif
