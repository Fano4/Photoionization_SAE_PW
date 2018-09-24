//
//  computation.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 20/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#include "algebra.hpp"


double determinant(double *A,int dim)
{
   //std::cout<<"probe2 dimension="<<dim<<std::endl;
    double det_val(1);
    short int sign(1);
    double *B=new double[dim*dim];
    int *ipiv=new int[dim];
    //double B[(dim-1)*(dim-1)];
    int n(0);
    for (int i=0; i!=dim; i++)
    {
        for (int j=0; j!=dim; j++)
        {
            B[i*dim+j]=A[i*dim+j];
        }
    }

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR,dim,dim,B,dim,ipiv);

/*    for (int i=0; i!=dim; i++)
    {
        for (int j=0; j!=dim; j++)
        {

            std::cout<<B[i*dim+j]<<"    ";
        }std::cout<<std::endl;
    }std::cout<<std::endl;*/
    for(int i=0;i!=dim;i++)
    {
       if(i+1!=ipiv[i])
          sign*=-1;

       det_val*=B[i*dim+i];

  //     std::cout<<i+1<<","<<ipiv[i]<<std::endl;
    }
//    std::cout<<"sign is "<<sign<<" ; det value is "<<det_val*sign<<std::endl;
     /*
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
                        //std::cout<<B[j*(dim-1)+n]<<"   ";
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
       // std::cout<<"VALUE = "<<A[0]<<std::endl;
        return A[0];
    }
    */
    delete [] ipiv;
    delete [] B;
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

double vector_prod(double vector1[],double vector2[],int gsize)
{
    double sum(0.);
//#pragma omp parallel for
        for (int j=0; j<gsize; j++)
        {
            sum+=vector1[j]*vector2[j];
        }

    return sum;
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
/*bool two_cubes_moment(double *cube1,double *cube2,double *moment,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
   double x(0);
   double y(0);
   double z(0);
   double dx((xmax-xmin)/nx);
   double dy((ymax-ymin)/ny);
   double dz((zmax-zmin)/nz);

   for(int i=0;i!=3;i++)
   {
      moment[i]=0;
   }
   for(int i=0;i!=nx;i++)
   {
      x=xmin+i*dx;
      for(int j=0;j!=ny;j++)
      {
         y=ymin+j*dy;
         for(int k=0;k!=nz;k++)
         {
           z=zmin+k*dz;
           moment[0]+=x*cube1[i*ny*nz+j*nz+k]*cube2[i*ny*nz+j*nz+k]*dx*dy*dz;
           moment[1]+=y*cube1[i*ny*nz+j*nz+k]*cube2[i*ny*nz+j*nz+k]*dx*dy*dz;
           moment[2]+=z*cube1[i*ny*nz+j*nz+k]*cube2[i*ny*nz+j*nz+k]*dx*dy*dz;
         }
      }
   }
   return 1;
}*/

double cube_dot_product(double *cube1,double *cube2,int nx,int ny, int nz,double dx,double dy,double dz,int angle_vec_size,double *output)
{
   double sum(0);
   int num(nx*ny*nz);
   int inc(1);
   const MKL_INT lda(num);
   int i(0);
   int j(0);

   //cblas_dgemv (CblasRowMajor, CblasNoTrans, angle_vec_size, num, dx*dy*dz ,cube1,lda, cube2, 1, 0,  output, 1);
  // #pragma omp parallel for private(i,j,sum) shared(angle_vec_size,nx,ny,nz,cube1,cube2,output,dx,dy,dz)
 /*  for(int i=0;i<angle_vec_size;i++)
   {
      output[i]=0;
      sum=0;
      for(int j=0;j<nx*ny*nz;j++)
      {
         sum+=cube1[i*nx*ny*nz+j]*cube2[j];
      }
      output[i]=sum*dx*dy*dz;
      //      sum+=cube1[i]*cube2[i]*dx*dy*dz;
      //output[i]= ddot(&num,&cube1[i*nx*ny*nz],&inc,cube2,&inc)*dx*dy*dz;
      // std::cout<<cube1[i*ny*nz+j*nz+k]<<" ; "<<cube2[i*ny*nz+j*nz+k]<<" ; "<<sum<<std::endl;
   }*/
   //std::cout<<sum<<" is the result of vector prod"<<std::endl;
   //return 0;
   *output= ddot(&num,cube1,&inc,cube2,&inc)*dx*dy*dz;
   return *output;
}
