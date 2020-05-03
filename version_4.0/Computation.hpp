#ifndef Computation_hpp
#define Computation_hpp

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include "gaunt_coeff_val.hpp"



bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff);
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
double MO_value( int mo_index, double r, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size,int **angular_mom_numbers);
double AO_value(int ao_index,double r, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int** angular_mom_numbers);
double contraction_value( double r, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type,int* angular_mom_numbers);
bool build_ao_s(double* S,int *nucl_basis_func,int *contraction_number,double **nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int basis_size);
void build_transition_density_matrix(int n_states_neut,int n_closed,int n_occ,int ci_size_neut,int n_elec_neut,double **ci_vec_neut,double **tran_den_mat_mo);
double build_reduced_determinant( int ai,int aj,int n_elec,int n_closed,int n_occ,double* mo_vector_1,double* mo_vector_2,double *spin_vector_1,double *spin_vector_2);
void compute_bessel_pice_mo(double*** pice_ortho_mo,double*** pice_ddx_mo,double*** pice_ddy_mo,double*** pice_ddz_mo,int jl_max,int n_occ,int basis_size,int nk,double kmax,double *MO_coeff_neutral,double **contraction_zeta,double **contraction_coeff,int * contraction_number,double** nucl_spher_pos,int *nucl_basis_func,int** angular_mom_numbers);


double associated_legendre(unsigned int l,int m,double x);
double associated_legendre_nonorm(unsigned int l,int m,double x);
double associated_legendre_der(unsigned int l,int m,double x);
//double legendre(unsigned int l,double x);
double legendre_der(unsigned int l,double x);
double rYlm (int l,int m,double thet,double phi);
double Dint(int l1,int l2,int l3,int m1,int m2,int m3);
double gaunt_formula(int l1,int l2,int l3,int m1,int m2,int m3);
void Jint_sort_indices(int* l1,int* l2,int* l3,int* m1,int* m2,int* m3);
double Jint_signflip_renormalize(int l1,int l2,int l3,int* m1,int* m2,int* m3);
double Jint_normalize(int l1,int l2,int l3,int m1,int m2,int m3);
double ALP_integral(int l,int m);
bool Jint_special_cases(int l1,int l2,int l3,int m1,int m2,int m3,double* result);
double prefactor_rYlm(int l,int m);
double azim_integ(int m1,int m2,int m3);
double polar_integ(int l1,int l2,int l3,int m1,int m2,int m3);
double I_p1_D_integral(int m1,int m2,int m3);
double I_m1_D_integral(int m1,int m2,int m3);
double I_p1_integral(int m1,int m2,int m3);
double I_m1_integral(int m1,int m2,int m3);
double J_int_D(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_p1_D(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1_D(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_p1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3);

#include "Legendre_functions.cpp"
#endif /*Computation_hpp*/
