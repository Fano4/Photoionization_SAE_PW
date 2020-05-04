#include <iostream>
#include <cmath>
#include <gsl/gsl_sf.h>
#include "prime.hpp"

double azim_integ(int m1,int m2,int m3)
{
   using namespace std;
   if(m1 ==0 && m2 ==0 && m3 ==0)
      return 2*acos(-1);
   else if(m1<0 && m2<0 && m3<0)
      return 0;
   else if(m1<0 && m2<0 && m3>=0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0)+bool(fabs(m1)-fabs(m2)+fabs(m3)==0)-bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else if(m1<0 && m2>=0 && m3<0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0)-bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else if(m1<0 && m2>=0 && m3>=0)
      return 0;
   else if(m1>=0 && m2<0 && m3<0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0)-bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0));
   else if(m1>=0 && m2<0 && m3>=0)
      return 0;
   else if(m1>=0 && m2>=0 && m3<0)
      return 0;

   else if(m1 >=0 && m2>=0 && m3 >=0)
      return 0.5*acos(-1)*(bool(fabs(m1)+fabs(m2)+fabs(m3)==0)+bool(fabs(m1)-fabs(m2)-fabs(m3)==0)+bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else 
   {
      std::cout<<"ERROR AZIM INTEGRAL TYPE CANNOT BE DEFINED "<<m1<<","<<m2<<","<<m3<<std::endl;
      exit(EXIT_SUCCESS);
   }
}
void Jint_sort_indices(int* l1,int* l2,int* l3,int* m1,int* m2,int* m3)
{
   //rearrange the array for getting l1<l2<l3
   int tmpl;
   int tmpm;
   if(*l1 > *l2)
   {
      tmpl=*l1;
      tmpm=*m1;
      *l1=*l2;
      *m1=*m2;
      *l2=tmpl;
      *m2=tmpm;
   }
   if(*l1 > *l3)
   {
      tmpl=*l1;
      tmpm=*m1;
      *l1=*l3;
      *m1=*m3;
      *l3=tmpl;
      *m3=tmpm;
   }
   if(*l2 > *l3)
   {
      tmpl=*l2;
      tmpm=*m2;
      *l2=*l3;
      *m2=*m3;
      *l3=tmpl;
      *m3=tmpm;
   }

   return;
}
double Jint_signflip_renormalize(int l1,int l2,int l3,int* m1,int* m2,int* m3)
{

   bool sgnm1( bool ( *m1 < 0 ) );
   bool sgnm2( bool ( *m2 < 0 ) );
   bool sgnm3( bool ( *m3 < 0 ) );
   double temp;
   double val(1);

   //Flip the sign of m's if they are negative
   if(sgnm1)
      *m1=-*m1;

   if(sgnm2)
      *m2=-*m2;

   if(sgnm3)
      *m3=-*m3;

   temp=1;
   if(sgnm1)
   {
      for (int tt=l1-*m1+1;tt!=l1+*m1+1;tt++)
         temp*=double(tt);
   }
   val/=temp;

   temp=1;
   if(sgnm2)
   {
         for (int tt=l2-*m2+1;tt!=l2+*m2+1;tt++)
            temp*=double(tt);
   }
   val/=temp;

   temp=1;
   if(sgnm3)
   {
      for (int tt=l3-*m3+1;tt!=l3+*m3+1;tt++)
          temp*=double(tt);
   }
   val/=temp;

   return pow(-1,sgnm1**m1+sgnm2**m2+sgnm3**m3)*val;

}
double Jint_normalize(int l1,int l2,int l3,int m1,int m2,int m3)
{
   double temp;
   double val(1);

      temp=1;
      for (int tt=l1-m1+1;tt!=l1+m1+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

      temp=1;
      for (int tt=l2-m2+1;tt!=l2+m2+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

      temp=1;
      for (int tt=l3-m3+1;tt!=l3+m3+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

      temp=val;
      return temp;
}
bool Jint_special_cases(int l1,int l2,int l3,int m1,int m2,int m3,double* result)
{
   int delta(m2+(m1-m2)*bool(m2>=m1));
   double temp;
   double prefactor(0);

   // Checking if any of the polynomials is zero
   
   if(l1 <0 || l2<0 || l3<0) //Negative degree should not occur
   {
      *result=0;
      return 1;
   }
   else if(m1>l1 || m2>l2 || m3>l3) // These ensure the degree of the polynomial is at least zero
   {
      *result=0;
      return 1;
   }

   else if( (l1+l2+l3) % 2 != (m1+m2+m3) % 2 ) // The integrand should be even,meaning that the sums of orders and degree have same parity
   {
      *result=0;
      return 1;
   }

   //If none of the polynomials is zero and the integrand is even,
   //
   // Sorting, sign flips and normalization
   
   Jint_sort_indices(&l1,&l2,&l3,&m1,&m2,&m3); // rearrange to get l1 < l2 < l3
   prefactor=Jint_signflip_renormalize(l1,l2,l3,&m1,&m2,&m3); // Flip the sign of negative m's and renormalize accordingly
   prefactor*=Jint_normalize(l1,l2,l3,m1,m2,m3); // Normalize the product of ALPs
   
   //Then check for special cases

   if( l1 == 0 ) //Overlap integral of two ALPs
   {
      temp=1;
      for (int tt = l2-m2+1;tt!=l2+m2+1;tt++)
      {
          temp*=double(tt);
      }
      *result = bool( l2 == l3) * 2 * temp / (2 * l2 + 1);
      return 1;
   }
   else if(m1+m2==m3)
   {
      *result = 2. * pow(-1,m3) * prefactor * gsl_sf_coupling_3j(2*l1,2*l2,2*l3,0,0,0)* gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,-2*m3);
      return 1;
   }

   else if(abs(m1+m2)==m3)
   {
      *result = 2.*pow(-1,delta-m1+m2)*prefactor* gsl_sf_coupling_3j(2*l1,2*l2,2*l3,0,0,0) * gsl_sf_coupling_3j(2*l1,2*l2,2*l3,-2*m1,2*m2,2*m1-2*m2);
      return 1;
   }
   else
      return 0;

}
double ALP_integral(int l,int m)
{

   if( (l+m) % 2 != 0)//Check for parity
      return 0;

   else if( l == 0 ) //Check for zero degree
      return 2;
   else
   {
      return (pow(double(-1),m)+pow(double(-1),l))
         *pow(2.,double(m-2))*double(m)
         *std::tgamma(double(l)/2.)
         *std::tgamma(double(l+m+1)/2.)
         /(std::tgamma(double(l-m)/2.+1)*std::tgamma(double(l+3)/2.));
   }
}
