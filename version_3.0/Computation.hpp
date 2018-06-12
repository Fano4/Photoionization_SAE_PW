#ifndef Computation_hpp
#define Computation_hpp

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "global_vars.hpp"
#include <complex>



bool Runge_kutta(double Rey0[], double Imy0[], double h, double time, int n_states, double *pot, double *dipole_x, double *dipole_y, double *dipole_z);
bool t_deriv(double* Repsi, double* Impsi, double* Redpsi, double* Imdpsi, double time, int n_states,double *pot, double *dipole_x, double *dipole_y, double *dipole_z);
bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff);

#endif /*Computation_hpp*/
