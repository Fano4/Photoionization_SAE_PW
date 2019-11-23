//
//  computation.hpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 20/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#ifndef computation_hpp
#define computation_hpp

#include <stdio.h>
#include <iostream>

void matrix_product(double *result_matrix,double *first_matrix,double *second_matrix,int dim1,int dim2,int dim3);
void transpose(double *input,double *transposed, int n_rows_input, int n_col_input);
double determinant(double *A,int dim);
//unsigned long long int factorial(int n,unsigned long long int* fact_memo);
bool kronecker_delta(int a, int b);
bool two_cubes_moment(double *cube1,double *cube2,double *moment,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
double cube_dot_product(double *cube1,double *cube2,int nx,int ny, int nz,double dx,double dy,double dz,int angle_vec_size,double *output);
double vector_prod(double vector1[],double vector2[],int gsize);
long double intplushalf_gamma(int n); //(Gamma(n+1/2))
long double gamma_int_or_half(double z);
double w3j(int l1,int l2,int l3,int m1,int m2,int m3);
double wdelta(int a,int b,int c);
double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3);
double j_l(int l,double z);//spherical bessel function of order l
double dj_ldz(int l,double z); //Derivative of the spherical bessel function of order l
int dfactorial(int n);
i//double ln_factorial(int n,double *memo);
void fact_prime_decomposer(int N, int* N_prime);

#endif /* computation_hpp */
