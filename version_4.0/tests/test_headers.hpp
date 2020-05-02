#ifndef test_file_hpp
#define test_file_hpp


//TEST_ALGEBRA.CPP

bool test_determinant();
bool test_prime_decomposer();
bool test_factorized_sum();
bool test_wigner3j();


//angular_int_aux_test.cpp

bool azim_integ_test();
bool gaunt_formula_test();
bool J_int_m2_test();

#include "../photoion_comp.hpp"

#include "test_utility.cpp"
#include "algebra_test.cpp"
#include "Legendre_functions_test.cpp"
#include "angular_int_aux_test.cpp"
#endif
