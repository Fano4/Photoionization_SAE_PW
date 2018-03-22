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
#include "global_vars.hpp"

void matrix_product(double *result_matrix,double *first_matrix,double *second_matrix,int dim1,int dim2,int dim3);
void transpose(double *input,double *transposed, int n_rows_input, int n_col_input);
double determinant(double *A,int dim);
long int factorial(int n);
bool kronecker_delta(int a, int b);
bool two_cubes_moment(double *cube1,double *cube2,double *moment,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
double cube_dot_product(double *cube1,double *cube2,int nx,int ny, int nz,double dx,double dy,double dz,int angle_vec_size,double *output);
double vector_prod(double vector1[],double vector2[],int gsize);

#endif /* computation_hpp */
