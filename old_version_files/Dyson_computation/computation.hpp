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
long int factorial(int n);
bool kronecker_delta(int a, int b);


#endif /* computation_hpp */
