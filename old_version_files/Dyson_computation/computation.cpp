//
//  computation.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 20/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#include "computation.hpp"


double determinant(double *A,int dim)
{
        double det_val(0.);
    double B[(dim-1)*(dim-1)];
    int n(0);

    /*for (int i=0; i!=dim; i++)
    {
        for (int j=0; j!=dim; j++)
        {

            std::cout<<A[i*dim+j]<<"    ";
        }std::cout<<std::endl;
    }std::cout<<std::endl;*/

    if (dim!=1)
    {
        for (int i=0; i!=dim; i++)
        {
            //std::cout<<"compute determinant  of "<<std::endl;
            for (int j=0; j!=dim-1; j++)
            {
                n=0;
                for (int k=0; k!=dim; k++)
                {
                    if(k!=i)
                    {
                        B[j*(dim-1)+n]=A[(j+1)*dim+k];
                        //std::cout/* <<"B"<<j<<" "<<n<<"="*/<<B[j*(dim-1)+n]<<"   ";
                        n++;
                    }
                }//std::cout<<std::endl;
            }//std::cout<<" ==== "<<std::endl;
            det_val+=pow(-1, i)*A[i]*determinant(B, dim-1);
        }
        //std::cout<<"VALUE = "<<det_val<<std::endl;
    }
    else
    {
        //std::cout<<"VALUE = "<<A[0]<<std::endl;
        return A[0];
    }
    

    
    
    
    return det_val;
}

void matrix_product(double *C,double *A,double *B,int dim1,int dim2,int dim3)
{
    double ntemp;
    
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim3; j++)
        {
            ntemp=0;
            
            for (int k=0; k!=dim2; k++)
            {
                ntemp+=A[i*dim2+k]*B[k*dim3+j];
                //std::cout<<A[i*dim2+k]<<"   *   "<<B[k*dim3+j]<<"   =   "<<A[i*dim2+k]*B[k*dim3+j]<<"  -> "<<ntemp<<std::endl;
            }
            C[i*dim3+j]=ntemp;
        }
    }
    
}
void transpose(double *A,double *B, int dim1, int dim2)
{
    //std::cout<<"   dim1= "<<dim1<<"   dim2= "<<dim2<<std::endl;
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim2; j++)
        {
            B[j*dim1+i]=A[i*dim2+j];
            //std::cout<<i*dim2+j<<"->"<<A[i*dim2+j]<<std::endl;
        }//std::cout<<std::endl;
    }
}

long int factorial(int n)
{
    if (n>1)
    return n*factorial(n-1);
    
    else
        return 1;

}

bool kronecker_delta(int a, int b)
{
    if (a==b)
    {
        return 1;
    }
    else
        return 0;
}
