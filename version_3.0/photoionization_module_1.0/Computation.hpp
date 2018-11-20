#ifndef Computation_hpp
#define Computation_hpp

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
//#include <gsl/gsl_sf_hyperg.h>


bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff);
std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> AO_FT(int ao_index,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers);
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);

std::complex<double> MO_Fourier_transform_grad( int mo_index,int comp, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int** angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);

std::complex<double> AO_FT_grad(int ao_index,int comp,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int *angular_mom_numbers);
std::complex<double> contraction_FT_grad_k( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> contraction_FT_grad_thet( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> contraction_FT_grad_phi( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
double MO_value( int mo_index, double r, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size,int **angular_mom_numbers);
double AO_value(int ao_index,double r, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int** angular_mom_numbers);
double contraction_value( double r, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type,int* angular_mom_numbers);
bool build_ao_s(double* S,int *nucl_basis_func,int *contraction_number,double **nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int basis_size);

double legendre(unsigned int l,double x);
double associated_legendre(unsigned int l,int m,double x);
double associated_legendre_nonorm(unsigned int l,int m,double x);
double legendre_der(unsigned int l,double x);
double associated_legendre_der(unsigned int l,int m,double x);
#endif /*Computation_hpp*/
