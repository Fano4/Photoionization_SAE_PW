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
unsigned long long int factorial(int n)
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
   else
   {
      if (n==1 || n==0)
         return 1;
      else
         return n*factorial(n-1);
   }
}

long double intplushalf_gamma(int n) //(Gamma(n+1/2))
{
   int fac1[MAX_FACTORIAL_PRIME];
   int fac2[MAX_FACTORIAL_PRIME];

   fact_prime_decomposer(2*n,fac1);
   fact_prime_decomposer(n,fac2);


   double temp=1;
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
//      std::cout<<"-- "<<PRIME[i]<<" ; "<<fac1[i]<<" - "<<fac2[i]<<std::endl;
      temp*=pow(PRIME[i],fac1[i]-fac2[i]);
   }
   return sqrt(acos(-1))*temp/(pow(4,n));

//   return sqrt(acos(-1))* exp(ln_factorial(2*n,lnfact_memo)-ln_factorial(n,lnfact_memo))/(pow(4,n));

//   return sqrt(acos(-1))*double(factorial(2*n,fact_memo))/(pow(4,n)*double(factorial(n,fact_memo)));
}

long double gamma_int_or_half(double z)
{
   int fac1[MAX_FACTORIAL_PRIME];

   if(ceil(z)==floor(z) && z > 0)
   {
      fact_prime_decomposer(int(z)-1,fac1);
      long double temp=1;
      for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      {
         temp*=pow(PRIME[i],fac1[i]);
      }
      return temp;
//      return exp(ln_factorial(int(z)-1,lnfact_memo));
   }
   else if(ceil(2*z) == floor(2*z) && z > 0)
      return intplushalf_gamma(int(z-0.5));
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
/*
double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3)
{
   //Apply selection rules
   if(m1+m2+m3 !=0)
      return 0;
   else if( l3 > fabs(l1+l2) || l3 < fabs(l1-l2))
      return 0;
   else if ((l1+l2+l3)%2 !=0)
      return 0;

   //Apply special cases

   else if(l1==l2 && l3==0 && m1 == -m2)
      return pow(-1,l1-m1)/(sqrt(2*l1+1));

   //initiate the prime representation arrays
   int wprefac_fac[MAX_FACTORIAL_PRIME]={0};
   int wdelta_fac[MAX_FACTORIAL_PRIME]={0};
   int Zw3j[MAX_FACTORIAL_PRIME]={0};
   int Zpw3j[MAX_FACTORIAL_PRIME]={0};

   double result(1);
   double sign(pow(-1,l1-l2-m3));

   //Fill the arrays
   w3j_prefac(l1,l2,l3,m1,m2,m3,wprefac_fac); //must be divived by two at the end
   wdelta(l1,l2,l3,wdelta_fac); //must be divived by two at the end
//   sign*=w3j(l1,l2,l3,m1,m2,m3,Zw3j,Zpw3j);//Sign must be accounted for here


   int temp1(1);
   int temp2(1);
   int temp3(1);
   std::cout<<"prefac : ";
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      std::cout<<wprefac_fac[i]<<",";
   }std::cout<<std::endl;
   std::cout<<"wdelta : ";
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      std::cout<<wdelta_fac[i]<<",";
   }std::cout<<std::endl;
   std::cout<<"Z : ";
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      std::cout<<Zw3j[i]<<",";
   }std::cout<<std::endl;
   std::cout<<"Zp : ";
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      std::cout<<Zpw3j[i]<<",";
   }std::cout<<std::endl;

   //Combine the arrays and compute the resulting value

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      result*=pow(double(PRIME[i]),double(wdelta_fac[i]+wprefac_fac[i])/2+double(Zpw3j[i] - Zw3j[i]));

   return sign*result;
   return 0;
}
bool w3j_prefac(int l1,int l2,int l3,int m1,int m2,int m3,int* fac)
{
   int fac1[MAX_FACTORIAL_PRIME]={0};
   int fac2[MAX_FACTORIAL_PRIME]={0};
   int fac3[MAX_FACTORIAL_PRIME]={0};
   int fac4[MAX_FACTORIAL_PRIME]={0};
   int fac5[MAX_FACTORIAL_PRIME]={0};
   int fac6[MAX_FACTORIAL_PRIME]={0};

   fact_prime_decomposer(l1+m1,fac1);
   fact_prime_decomposer(l1-m1,fac2);
   fact_prime_decomposer(l2+m2,fac3);
   fact_prime_decomposer(l2-m2,fac4);
   fact_prime_decomposer(l3+m3,fac5);
   fact_prime_decomposer(l3-m3,fac6);

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      fac[i]=fac1[i]+fac2[i]+fac3[i]+fac4[i]+fac5[i]+fac6[i]; //!!! fac must be divived by two at the end

   return 0;
}
bool wdelta(int a,int b,int c,int* fac)
{

   int fac1[MAX_FACTORIAL_PRIME]={0};
   int fac2[MAX_FACTORIAL_PRIME]={0};
   int fac3[MAX_FACTORIAL_PRIME]={0};
   int fac4[MAX_FACTORIAL_PRIME]={0};

   fact_prime_decomposer(a+b-c,fac1);
   fact_prime_decomposer(a-b+c,fac2);
   fact_prime_decomposer(-a+b+c,fac3);
   fact_prime_decomposer(a+b+c+1,fac4);

   double temp(1);

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      fac[i]=fac1[i]+fac2[i]+fac3[i]-fac4[i];
   }

   return 0;
}
int w3j(int l1,int l2,int l3,int m1,int m2,int m3,int* Z,int* Zp)
{
   //We compute the value by using prime number factorization.
   //
   //We compute the expression using the form A=(1/Z)*sum_t (-1)**T * Z_t
   //
   //Z= ∏_t z**t
   //
   //Z_t= (-1)**t ∏_(t'!=t) z**t'
   //

   int fac1[MAX_FACTORIAL_PRIME]={0};
   int fac2[MAX_FACTORIAL_PRIME]={0};
   int fac3[MAX_FACTORIAL_PRIME]={0};
   int fac4[MAX_FACTORIAL_PRIME]={0};
   int fac5[MAX_FACTORIAL_PRIME]={0};
   int fac6[MAX_FACTORIAL_PRIME]={0};

   double sum(0);
   double val(0);
   double val2(0);
   double temp(1);

   //First, determine the range of the t index
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

   //Compute only if there is a valid range
   if(tmin <= tmax)
   {

      //initiate the arrays 
      
      int** Zpt=new int* [(tmax-tmin+1)];
      for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
         Z[i]=0;
      int sign(1);
      int temp_Z[MAX_FACTORIAL_PRIME]={0};
      for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
         Zp[i]=0;
      
      // Compute the prefactor Z

      for(int t=tmin;t!=tmax+1;t++)
      {
         std::cout<<"probe "<<tmax-tmin+1<<" : "<<t-tmin<<std::endl;
         fact_prime_decomposer(t,fac1);
         fact_prime_decomposer(l3-l2+t+m1,fac2);
         fact_prime_decomposer(l3-l1+t-m2,fac3);
         fact_prime_decomposer(l2+l1-t-l3,fac4);
         fact_prime_decomposer(l1-t-m1,fac5);
         fact_prime_decomposer(l2+t+m2,fac6);

         for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
         {
            Z[i]+=fac1[i]+fac2[i]+fac3[i]+fac4[i]+fac5[i]+fac6[i];
         }
      }

      //compute the factors Z_t

      for(int t=tmin;t!=tmax+1;t++)
      {
         Zpt[t-tmin]=new int[MAX_FACTORIAL_PRIME];
         std::cout<<"probe "<<t-tmin<<std::endl;

         for(int tp=tmin;tp!=tmax+1;tp++)
         {
            if(tp==t)
               continue;
            else
            {
               fact_prime_decomposer(tp,fac1);
               fact_prime_decomposer(l3-l2+tp+m1,fac2);
               fact_prime_decomposer(l3-l1+tp-m2,fac3);
               fact_prime_decomposer(l2+l1-tp-l3,fac4);
               fact_prime_decomposer(l1-tp-m1,fac5);
               fact_prime_decomposer(l2-tp+m2,fac6);
               for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
               {
                  std::cout<<"probe "<<t-tmin<<","<<i<<std::endl;
                  Zpt[t-tmin,i]+=fac1[i]+fac2[i]+fac3[i]+fac4[i]+fac5[i]+fac6[i];
               }std::cout<<"loop done"<<std::endl;
            }
         }

         //Now, add each Zpt, with the correct sign (-1)**t


         //sum up the factors Z_t with correct sign
               //first, identify the sign of both term
               // sign of the main term Zp is sign_Zp ; The sign of the additional term is (-1)**t
               //If the sign is different, then call difference
               //If the sign is the same, call addition ; The sign does not change

               if( t % 2 == 0 && sign > 0 )  // + + 
                  factorized_sum(Zp,Zpt[t-tmin],temp_Z);

               else if ( t % 2 == 0 && sign < 0 ) // + -
                  sign=factorized_diff(Zpt[t-tmin],Zp,temp_Z);

               else if ( t % 2 != 0 && sign > 0 ) // - +
                  sign=factorized_diff(Zp,Zpt[t-tmin],temp_Z);

               else // - -
                  factorized_sum(Zp,Zpt[t-tmin],temp_Z);

               //Keep track of the overal sign!
            
         for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
            Zp[i]=temp_Z[i]; 
      }
      //Eventually return the overal sign
      //
      //
      delete [] Zpt;
      return sign;
   }
   else
   {
      return 0;
   }

}
*/
double j_l(int l,double z) //spherical bessel function of order l
{
   double test(1);
   double val(0);
   double valm1(0);
   double factor(1);
   int i(0);
   int fac1[MAX_FACTORIAL_PRIME];
   double temp(1);


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
        else
        {
            val=1;
            i=1; 
            while(test>=1e-15)
            {
               fact_prime_decomposer(i,fac1);
               temp=1;
               for(int ii=0;ii!=MAX_FACTORIAL_PRIME;ii++)
               {
                  temp*=pow(PRIME[ii],fac1[ii]);
               }
                valm1=val;
                factor=1;
                for(int k=1;k!=i+1;k++) factor*=(2*k+2*l+1);
                val+=pow(-1,i)*pow(z*z/2.,i)/(factor*temp);
                i++;
                test=fabs((val-valm1));
            }

            if(isnan(val))
               std::cout<<" ERROR ! BESSEL FUNCTION IS NAN"<<std::endl;

            return val*pow(z,l)/dfactorial(2*l+1);
        }
   }
}

double dj_ldz(int l,double z) //Derivative of the spherical bessel function of order l
{
   if(z==0 || l == 0)
      return 0;
   else
   {
      if(isnan(l*j_l(l-1,z)-(l+1)*j_l(l+1,z))/(2*l+1))
          std::cout<<" ERROR ! BESSEL DERIVATIVE FUNCTION IS NAN"<<std::endl;
      return (l*j_l(l-1,z)-(l+1)*j_l(l+1,z))/(2*l+1);
   }
}
int dfactorial(int n)
{
   if(n<=1)
      return 1;
   else
      return n*dfactorial(n-2);
}
void prime_decomposer(int N, int* N_prime)
{
   //This routine factorizes an integer number into prime numbers.
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      N_prime[i]=0;

   if(N==0)
   {
      std::cout<<"WARNING ! TRYING TO DECOMPOSE ZERO IN PRIME NUMBERS"<<std::endl;
      return;
   }
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      while( N%PRIME[i] == 0 )
      {
         N/=PRIME[i];
         N_prime[i]++;
      }
   }
   if(N!=1)
   {
      std::cout<<" INCOMPLETE FACTORIZATION. Remaining factor : "<<N<<std::endl;
   }
   return;
}
void fact_prime_decomposer(int N, int* N_prime)
{
   //This routine factorizes the factorial of an integer number into prime numbers. It uses a global and constant array with maximum integer MAX_N_FACTORIAL
   if(N>MAX_N_FACTORIAL)
   {
      std::cout<<"ERROR ! FACTORIAL ARGUMENT LARGER THAN MAX AUTHORIZED VALUE ! N = "<<N<<std::endl;
      exit(EXIT_SUCCESS);
   }
   else
   {
      for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      {
         N_prime[i]=PRIME_DECOMPOSED_FAC[N][i];
      }
   } 
   return;
}
bool factorized_sum(int* x1,int* x2,int* out)
{
   //!!!!APPLIES ONLY TO FACTORIALS OF NUMBERS < MAX_N_FACTORIAL
   //We want to compute the factorized representation of a sum of two integers
   //1. we factorize the sum by comparing the vectors
   //2. We explicitly compute the remaining sum
   //3. we factorize the remaining sum and put everythingout
   
   int remain1(1);
   int remain2(1);
   int remain[MAX_FACTORIAL_PRIME];

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      out[i]=0;
      while(x1[i] != 0 && x2[i]!= 0)
      {
         out[i]++;
         x1[i]--;
         x2[i]--;
      }
      remain1*=pow(PRIME[i],x1[i]);
      remain2*=pow(PRIME[i],x2[i]);
   }
   prime_decomposer(remain1+remain2,remain);
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      out[i]+=remain[i];
   }
   return 0;
}
int factorized_diff(int* x1,int* x2,int* out)
{
   //!!!!APPLIES ONLY TO FACTORIALS OF NUMBERS < MAX_N_FACTORIAL
   //We want to compute the factorized representation of a difference of two integers x1-x2
   //1. we factorize the sum by comparing the vectors
   //2. We explicitly compute the remaining difference
   //3. we factorize the remaining difference and keep track of the sign
   //4. We put the resulting integer in the out array and the sign as a return
   
   int sign(1);
   int remain1(1);
   int remain2(1);
   int remain[MAX_FACTORIAL_PRIME];

   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      out[i]=0;
      while(x1[i] != 0 && x2[i]!= 0)
      {
         out[i]++;
         x1[i]--;
         x2[i]--;
      }
      remain1*=pow(PRIME[i],x1[i]);
      remain2*=pow(PRIME[i],x2[i]);
   }
   //check the sign of the difference and factorize the difference
   if(remain1 == remain2)
      return 0;
   else if(remain1>remain2)
      prime_decomposer(remain1-remain2,remain);
   else
   {
      prime_decomposer(remain2-remain1,remain);
      sign=-1;
   }
   for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
   {
      out[i]+=remain[i];
   }
   return sign;
}
