#ifndef Computation_hpp
#define Computation_hpp

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dawson.h>


bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff);
std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> AO_FT(int ao_index,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type);
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type);
std::complex<double> MO_Fourier_transform_grad( int mo_index,int comp, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> AO_FT_grad(int ao_index,int comp,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type);
std::complex<double> contraction_FT_grad_k( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type);
std::complex<double> contraction_FT_grad_thet( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type);
std::complex<double> contraction_FT_grad_phi( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type);

#endif /*Computation_hpp*/
