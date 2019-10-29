#ifndef Computation_hpp
#define Computation_hpp

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>



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
void build_transition_density_matrix(int n_states_neut,int n_closed,int n_occ,int ci_size_neut,int n_elec_neut,double **ci_vec_neut,double **tran_den_mat_mo);
double build_reduced_determinant( int ai,int aj,int n_elec,int n_closed,int n_occ,double* mo_vector_1,double* mo_vector_2,double *spin_vector_1,double *spin_vector_2);


double associated_legendre(unsigned int l,int m,double x);
double associated_legendre_nonorm(unsigned int l,int m,double x);
double associated_legendre_der(unsigned int l,int m,double x);
double legendre(unsigned int l,double x);
double legendre_der(unsigned int l,double x);
std::complex<double> bessel_MO_overlap( int mo_index, double k, int jl, int jml,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> bessel_AO_overlap(int ao_index,double k, int jl, int jml,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
double rYlm (int l,int m,double thet,double phi);
double Dint(int l1,int l2,int l3,int m1,int m2,int m3);
double gaunt_formula(int l1,int l2,int l3,int m1,int m2,int m3);
double prefactor_rYlm(int l,int m);
double azim_integ(int m1,int m2,int m3);
std::complex<double> bessel_azimder_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_polarder_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_radialder_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_ddx_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_ddx_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_ddy_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_ddz_contraction_overlap( double k, int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_AO_ddx_overlap(int ao_index,double k, int jl, int jml,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_AO_ddy_overlap(int ao_index,double k, int jl, int jml,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_AO_ddz_overlap(int ao_index,double k, int jl, int jml,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers);
std::complex<double> bessel_MO_ddx_overlap( int mo_index, double k, int jl, int jml,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> bessel_MO_ddy_overlap( int mo_index, double k, int jl, int jml,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);
std::complex<double> bessel_MO_ddz_overlap( int mo_index, double k, int jl, int jml,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size);

#include "Legendre_functions.cpp"
#endif /*Computation_hpp*/
