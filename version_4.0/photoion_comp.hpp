#ifndef FRONUC_DYN_M_H
#define FRONUC_DYN_M_H
#define MKL_Complex16 std::complex<double>
#define MAX_N_FACTORIAL 20
#define MAX_LN_FACTORIAL 50
#define MAX_FACTORIAL_PRIME 25

const int PRIME[MAX_FACTORIAL_PRIME]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

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
#include "files_reader.hpp"
#include "dyson_cube_writer.hpp"
#include "algebra.hpp"
#include "Computation.hpp"
#include "hf5_photoion.hpp"

#include "files_reader.cpp"
#include "dyson_cube_writer.cpp"
#include "algebra.cpp"
#include "angular_int_aux.cpp"
#include "Computation.cpp"
#include "hf5_photoion.cpp"
#include "test_file.cpp"

int spherical_harmonics_translator(std::string basis_func_type,bool component);

#endif
