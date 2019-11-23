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
    double det_val(1);
    short int sign(1);
    double *B=new double[dim*dim];
    int *ipiv=new int[dim];
    int n(0);
    for (int i=0; i!=dim; i++)
    {
        for (int j=0; j!=dim; j++)
        {
            B[i*dim+j]=A[i*dim+j];
        }
    }

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR,dim,dim,B,dim,ipiv);

    for(int i=0;i!=dim;i++)
    {
       if(i+1!=ipiv[i])
          sign*=-1;

       det_val*=B[i*dim+i];
    }
    delete [] ipiv;
    delete [] B;
    return sign*det_val;
}

void matrix_product(double *C,double *A,double *B,int dim1,int dim2,int dim3)
{
   //C=A*B
    double ntemp;
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim3; j++)
        {
            ntemp=0;
            
            for (int k=0; k!=dim2; k++)
            {
                ntemp+=A[i*dim2+k]*B[k*dim3+j];
            }
            C[i*dim3+j]=ntemp;
        }
    }
    
}
void transpose(double *A,double *B, int dim1, int dim2)
{
   //B=trans(A)
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim2; j++)
        {
            B[j*dim1+i]=A[i*dim2+j];
        }
    }
}
/*
double ln_factorial(int n,double *memo)
{
   if(n>MAX_LN_FACTORIAL)
   {
      std::cout<<"WARNING LARGE FACTORIAL ARGUMENT : N ="<<n<<std::endl<<"EXIT"<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else if(n<0)
   {
      std::cout<<"FATAL ERROR! NEGATIVE ARGUMENT IN FACTORIAL"<<std::endl<<"N = "<<n<<std::endl<<"EXIT"<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else if(!(memo[n]==0))
   {
      return memo[n];
   }
   else
   {
      if(n==0)
         memo[n]=0;
      else
         memo[n]=log(n)+ln_factorial(n-1,memo);
      return memo[n];
   }
}
unsigned long long int factorial(int n,unsigned long long int* memo)
{
   if(n>MAX_N_FACTORIAL)
   {
      std::cout<<"WARNING LARGE FACTORIAL ARGUMENT : N ="<<n<<std::endl<<"EXIT"<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else if(n<0)
   {
      std::cout<<"FATAL ERROR! NEGATIVE ARGUMENT IN FACTORIAL"<<std::endl<<"N = "<<n<<std::endl<<"EXIT"<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else if(!(memo[n]==0))
   {

      return memo[n];
   }
   else
   {
     // std::cout<<"Computing factorial of "<<n<<std::endl;
      if(n==0)
         memo[n]=1;
      else
         memo[n]=n*factorial(n-1,memo);
      return memo[n];
   }
}
*/
long double intplushalf_gamma(int n,double* lnfact_memo) //(Gamma(n+1/2))
{
   return sqrt(acos(-1))* exp(ln_factorial(2*n,lnfact_memo)-ln_factorial(n,lnfact_memo))/(pow(4,n));
//   return sqrt(acos(-1))*double(factorial(2*n,fact_memo))/(pow(4,n)*double(factorial(n,fact_memo)));
}

long double gamma_int_or_half(double z,double* lnfact_memo)
{
   if(ceil(z)==floor(z) && z > 0)
      return exp(ln_factorial(int(z)-1,lnfact_memo));
   else if(ceil(2*z) == floor(2*z) && z > 0)
      return intplushalf_gamma(int(z-0.5),lnfact_memo);
   else
   {
      std::cout<<"ERROR ! NON HALF INTEGER OR NON INTEGER OR NON POSITIVE ARGUMENT IN GAMMA_INT_OR_HALF FUCNTION"<<std::endl;
      std::cout<<" Z = "<<z<<std::endl;
      exit(EXIT_SUCCESS);
   }
}

double vector_prod(double vector1[],double vector2[],int gsize)
{
    double sum(0.);
#pragma omp parallel for
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
double cube_dot_product(double *cube1,double *cube2,int nx,int ny, int nz,double dx,double dy,double dz,int angle_vec_size,double *output)
{
   double sum(0);
   int num(nx*ny*nz);
   int inc(1);
   const MKL_INT lda(num);

   #pragma omp parallel for 
   for(int i=0;i<angle_vec_size;i++)
   {
      output[i]=0;
      sum=0;
      for(int j=0;j<nx*ny*nz;j++)
      {
         sum+=cube1[i*nx*ny*nz+j]*cube2[j];
      }
      output[i]=sum*dx*dy*dz;
   }
   return 0;
}

double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   double val(pow(-1,l1-l2-m3)*wdelta(l1,l2,l3,lnfact_memo)*w3j(l1,l2,l3,m1,m2,m3,lnfact_memo));

/*   if(val == 0 )
   {
      std::cout<<"Inside Wigner"<<std::endl;
      std::cout<<wdelta(l1,l2,l3,fact_memo)<<std::endl;
      std::cout<<w3j(l1,l2,l3,m1,m2,m3,fact_memo)<<std::endl;
   }*/

   if(isnan(val))
       std::cout<<" ERROR ! WIGNER3J FUNCTION IS NAN"<<std::endl;
   return val;
}
double wdelta(int a,int b,int c,double* lnfact_memo)
{

   //in order to increase the precision of the evaluation of the facotorial, we will decompose all factorials into it's prime numbers. We wil then simplify the product of factorials by computing the sum of the powers of the prime numbers before to compute the final value.
   //
   //1. What is the largest of the numbers for which we compute the exponential and how many prime numbers are contained between 1 and this number? This gives the size of the array.
   //2. we build as many arrays as there are factorials to evaluate and we decompose these facotrials.
   //3. We sum the powers of the arrays to obtain the final result. Then we compute the actual result.

   
   double val(exp(0.5*(ln_factorial(a+b-c,lnfact_memo)+ln_factorial(a-b+c,lnfact_memo)+ln_factorial(-a+b+c,lnfact_memo)-ln_factorial(a+b+c+1,lnfact_memo))));
   //double val( sqrt(double(factorial(a+b-c,fact_memo))*double(factorial(a-b+c,fact_memo))*double(factorial(-a+b+c,fact_memo))/double(factorial(a+b+c+1,fact_memo))));
   if(isnan(val))
       std::cout<<" ERROR ! WDelta FUNCTION IS NAN"<<std::endl;
   return val;
}
double w3j(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   double sum(0);
   double val(0);
   double val2(0);
   double temp(1);

   int tmin(0);
   int tmax(1000);

   if(l2-l3-m1>tmin)
      tmin=l2-l3-m1;
   if(l1-l3+m2>tmin)
      tmin=l1-l3+m2;
   
   if(l1-m1<tmax)
      tmax=l1-m1;
   if(l2+m2<tmax)
      tmax=l2+m2;
   if(l2+l1-l3<tmax)
      tmax=l2+l1-l3;

//   std::cout<<"numbers"<<std::endl<<l1<<","<<l2<<","<<l3<<","<<m1<<","<<m2<<","<<m3<<std::endl;
//   std::cout<<tmin<<","<<tmax<<std::endl;
   if(tmin <= tmax)
   {
      for(int t=tmin;t!=tmax+1;t++)
      {
         sum+=pow(-1.,t) * exp(-(ln_factorial(t,lnfact_memo) + ln_factorial(l3-l2+t+m1,lnfact_memo) + ln_factorial(l3-l1+t-m2,lnfact_memo) + ln_factorial(l2+l1-t-l3,lnfact_memo) + ln_factorial(l1-t-m1,lnfact_memo) + ln_factorial(l2-t+m2,lnfact_memo)));
//         sum+=pow(-1.,t)/(double(factorial(t,fact_memo))*double(factorial(l3-l2+t+m1,fact_memo))*double(factorial(l3-l1+t-m2,fact_memo))*double(factorial(l2+l1-t-l3,fact_memo))*double(factorial(l1-t-m1,fact_memo))*double(factorial(l2-t+m2,fact_memo)));
         if(isnan(sum))
            std::cout<<t<<","<<l3-l2+t+m1<<","<<(l3-l1+t-m2)<<","<<l2+l1-t-l3<<","<<l1-t-m1<<","<<l2-t+m2<<std::endl;
      }
   }
   else
   {
      return 0;
   }

//   std::cout<<"====="<<sum<<std::endl;
      temp=1;
      val2=1;
      for (int tt=l3-m3+1;tt!=l3+m3+1;tt++)
      {
         temp*=tt;
      }
      val2*=sqrt(temp);
   val=exp(0.5*(ln_factorial(l1+m1,lnfact_memo) + ln_factorial(l1-m1,lnfact_memo) + ln_factorial(l2+m2,lnfact_memo) + ln_factorial(l2-m2,lnfact_memo) + ln_factorial(l3+m3,lnfact_memo) + ln_factorial(l3-m3,lnfact_memo)))*sum;
//   val=sqrt(double(factorial(l1+m1,fact_memo)))*sqrt(double(factorial(l1-m1,fact_memo)))*sqrt(double(factorial(l2+m2,fact_memo)))*sqrt(double(factorial(l2-m2,fact_memo)))*sqrt(double(factorial(l3+m3,fact_memo)))*sqrt(double(factorial(l3-m3,fact_memo)))*sum;
   if(isnan(val))
   {
       std::cout<<" ERROR ! W3j FUNCTION IS NAN"<<std::endl<<l1<<","<<l2<<","<<l3<<";"<<m1<<","<<m2<<","<<m3<<std::endl;
       std::cout<<sum<<std::endl;
//       std::cout<<factorial(l1+m1,fact_memo)<<","<<factorial(l1-m1,fact_memo)<<","<<factorial(l2+m2,fact_memo)<<","<<factorial(l2-m2,fact_memo)<<","<<factorial(l3+m3,fact_memo)<<","<<factorial(l3-m3,fact_memo)<<std::endl;
       exit(EXIT_SUCCESS);
   }
   return val;
}
double j_l(int l,double z,double* lnfact_memo) //spherical bessel function of order l
{
   double test(1);
   double val(0);
   double valm1(0);
   double factor(1);
   int i(0);

   if(z==0)
   {
      if(l==0)
         return 1;
      else
         return 0;
   }
   else
   {
        if(l<-1)
            return 0;
        else if(l == -1)
           return cos(z)/z;
//        else if(l==0)
//            return 1;
        else
        {
            val=1;
            i=1; 
            while(test>=1e-15)
            {
                valm1=val;
                factor=1;
                for(int k=1;k!=i+1;k++) factor*=(2*k+2*l+1);
                val+=pow(-1,i)*pow(z*z/2,i)/(factor*exp(ln_factorial(i,lnfact_memo)));
                i++;
                test=fabs((val-valm1));
            }

            if(isnan(val))
               std::cout<<" ERROR ! BESSEL FUNCTION IS NAN"<<std::endl;

            return val*pow(z,l)/dfactorial(2*l+1);
        }
   }
}

double dj_ldz(int l,double z,double* lnfact_memo) //Derivative of the spherical bessel function of order l
{
   if(z==0 || l == 0)
      return 0;
   else
   {
      if(isnan(l*j_l(l-1,z,lnfact_memo)-(l+1)*j_l(l+1,z,lnfact_memo))/(2*l+1))
          std::cout<<" ERROR ! BESSEL DERIVATIVE FUNCTION IS NAN"<<std::endl;
      return (l*j_l(l-1,z,lnfact_memo)-(l+1)*j_l(l+1,z,lnfact_memo))/(2*l+1);
   }
}
int dfactorial(int n)
{
   if(n<=1)
      return 1;
   else
      return n*dfactorial(n-2);
}
void fact_prime_decomposer(int N, int* N_prime)
{
   int prime[25]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
   int m=0;

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      N_prime[i]=0;
      for(int n=2;n!=N+1;n++)
      {
         m=n;
         while(m%prime[i]==0)
         {
            N_prime[i]++;
            m/=N_prime[i];
         }
      }
   }

}
