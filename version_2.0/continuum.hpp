#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "global_vars.hpp"
#include "dyson_cube_writer.hpp"
#include "algebra.hpp"

bool pw_builder(double *Reinput,double *Iminput,double k,double *theta,double *phi,int angle_vec_size,double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz,int n_occ,double **neut_mo_cube_array);
bool pw_orthogonalizer(double *Reinput,double *Iminput,int angle_vec_size,int nx,int ny,int nz,double dx,double dy,double dz,int n_occ,double **neut_mo_cube_array);
bool projector_lz(double *Reinput,double *Iminput,int angle_vec_size,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,int lz);
