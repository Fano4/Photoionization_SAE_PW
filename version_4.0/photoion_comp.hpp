#ifndef FRONUC_DYN_M_H
#define FRONUC_DYN_M_H
#define MKL_Complex16 std::complex<double>
#define MAX_N_FACTORIAL 100
#define MAX_LN_FACTORIAL 50
#define MAX_FACTORIAL_PRIME 25

const int PRIME[MAX_FACTORIAL_PRIME]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
const int PRIME_DECOMPOSED_FAC[MAX_N_FACTORIAL+1][MAX_FACTORIAL_PRIME]={
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{3,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{4,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{4,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{7,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{7,4,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{8,4,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{8,4,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{10,5,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{10,5,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{11,5,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{11,6,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{15,6,3,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{15,6,3,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{16,8,3,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{16,8,3,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{18,8,4,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{18,9,4,3,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{19,9,4,3,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{19,9,4,3,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{22,10,4,3,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{22,10,6,3,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{23,10,6,3,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{23,13,6,3,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{25,13,6,4,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{25,13,6,4,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{26,14,7,4,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{26,14,7,4,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{31,14,7,4,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{31,15,7,4,3,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{32,15,7,4,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{32,15,8,5,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{34,17,8,5,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{34,17,8,5,3,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
{35,17,8,5,3,2,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
{35,18,8,5,3,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
{38,18,9,5,3,3,2,2,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
{38,18,9,5,3,3,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
{39,19,9,6,3,3,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
{39,19,9,6,3,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0},
{41,19,9,6,4,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0},
{41,21,10,6,4,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0},
{42,21,10,6,4,3,2,2,2,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0},
{42,21,10,6,4,3,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{46,22,10,6,4,3,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{46,22,10,8,4,3,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{47,22,12,8,4,3,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{47,23,12,8,4,3,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{49,23,12,8,4,4,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0},
{49,23,12,8,4,4,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{50,26,12,8,4,4,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{50,26,13,8,5,4,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{53,26,13,9,5,4,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{53,27,13,9,5,4,3,3,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{54,27,13,9,5,4,3,3,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0},
{54,27,13,9,5,4,3,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0},
{56,28,14,9,5,4,3,3,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0},
{56,28,14,9,5,4,3,3,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{57,28,14,9,5,4,3,3,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{57,30,14,10,5,4,3,3,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{63,30,14,10,5,4,3,3,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{63,30,15,10,5,5,3,3,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{64,31,15,10,6,5,3,3,2,2,2,1,1,1,1,1,1,1,0,0,0,0,0,0,0},
{64,31,15,10,6,5,3,3,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0},
{66,31,15,10,6,5,4,3,2,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0},
{66,32,15,10,6,5,4,3,3,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0},
{67,32,16,11,6,5,4,3,3,2,2,1,1,1,1,1,1,1,1,0,0,0,0,0,0},
{67,32,16,11,6,5,4,3,3,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0},
{70,34,16,11,6,5,4,3,3,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0,0},
{70,34,16,11,6,5,4,3,3,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0,0},
{71,34,16,11,6,5,4,3,3,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0},
{71,35,18,11,6,5,4,3,3,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0},
{73,35,18,11,6,5,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0},
{73,35,18,12,7,5,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0},
{74,36,18,12,7,6,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0,0},
{74,36,18,12,7,6,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0},
{78,36,19,12,7,6,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0},
{78,40,19,12,7,6,4,4,3,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0,0},
{79,40,19,12,7,6,4,4,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0,0},
{79,40,19,12,7,6,4,4,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0},
{81,41,19,13,7,6,4,4,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0},
{81,41,20,13,7,6,5,4,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0,0},
{82,41,20,13,7,6,5,4,3,2,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0},
{82,42,20,13,7,6,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0},
{85,42,20,13,8,6,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0,0},
{85,42,20,13,8,6,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0},
{86,44,21,13,8,6,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0},
{86,44,21,14,8,7,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0},
{88,44,21,14,8,7,5,4,4,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1,0},
{88,45,21,14,8,7,5,4,4,3,3,2,2,2,1,1,1,1,1,1,1,1,1,1,0},
{89,45,21,14,8,7,5,4,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0},
{89,45,22,14,8,7,5,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0},
{94,46,22,14,8,7,5,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,0},
{94,46,22,14,8,7,5,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1},
{95,46,22,16,8,7,5,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1},
{95,48,22,16,9,7,5,5,4,3,3,2,2,2,2,1,1,1,1,1,1,1,1,1,1}};


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
//#include "dyson_cube_writer.hpp"
#include "algebra.hpp"
#include "Computation.hpp"
#include "hf5_photoion.hpp"

#include "files_reader.cpp"
//#include "dyson_cube_writer.cpp"
#include "algebra.cpp"
#include "angular_int_aux.cpp"
#include "Computation.cpp"
#include "hf5_photoion.cpp"
#include "test_file.cpp"

int spherical_harmonics_translator(std::string basis_func_type,bool component);

#endif