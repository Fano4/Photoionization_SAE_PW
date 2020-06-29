//testing routines for angular_int_aux.cpp
//
//
//

#include <iostream>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include "mathfunctions.h"
#include "prime.hpp"

//////////////////////////////////////////
//
//////////////////////////////////////////
bool two_azim_integ_test()
{
   std::cout<<"Testing two_azim_integ"<<std::endl;
   int m1,m2;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);

   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral gives good values
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=1000;i++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*phi(m1,x[i])*phi(m2,x[i]);
       
       //generate function result
       val=two_azim_integ(m1,m2);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//////////////////////////////////////////
bool three_azim_integ_test()
{
   std::cout<<"Testing three_azim_integ"<<std::endl;
   int m1,m2,m3;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);
   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral give good values
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=1000;i++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);
       m3=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*phi(m1,x[i])*phi(m2,x[i])*phi(m3,x[i]);
       
       //generate function result
       val=three_azim_integ(m1,m2,m3);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<","<<m3<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_Jint_sort_indices()
{
   bool test1(1);
   bool test2(1);
   std::cout<<"Testing Jint_sort_indices"<<std::endl;
   int m1,m2,m3,l1,l2,l3;

   //test 1 : Get the right order
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      //generate a random set of values
       l1=(rand() % 50);
       l2=(rand() % 50);
       l3=(rand() % 50);
       m1=( -l1 + rand() % (2*l1+2));
       m2=( -l2 + rand() % (2*l3+2));
       m3=( -l3 + rand() % (2*l2+2));

       //sort the values 
       Jint_sort_indices(&l1,&l2,&l3,&m1,&m2,&m3);

       //Check if the order is right
       if(l1>l2 || l1>l3 || l2>l3)
          test1*=0;
       else
          test1*=1;
   }

   //test 2 : Keep track of the values
   std::cout<<"2..."<<std::endl;
   {
      int cl1,cl2,cl3,cm1,cm2,cm3;
   for(int i=0;i!=25;i++)
   {
      //generate a random set of values
       l1=(rand() % 50);
       l2=(rand() % 50);
       l3=(rand() % 50);
       m1=( -l1 + rand() % (2*l1+2));
       m2=( -l2 + rand() % (2*l3+2));
       m3=( -l3 + rand() % (2*l2+2));

       //Copy the values
       cl1=l1;
       cl2=l2;
       cl3=l3;
       cm1=m1;
       cm2=m2;
       cm3=m3;

       //sort the values 
       Jint_sort_indices(&l1,&l2,&l3,&m1,&m2,&m3);

       //Check if anyone has been lost on the way
       if( (cl1 == l1 || cl1==l2 || cl1==l3) && (cl2 == l1 || cl2==l2 || cl2==l3) && ( cl3 == l1 || cl3==l2 || cl3==l3)  
             && (cm1 == m1 || cm1==m2 || cm1==m3) && (cm2 == m1 || cm2==m2 || cm2==m3) && (cm3 == m1 || cm3==m2 || cm3==m3) )
          test2*=1;
       else
          test2*=0;
   }
   }
   if(test1 && test2)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      if(!test2)
         std::cout<<"Error 2...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_Jint_signflip_renormalize()
{
   bool test1(1);
   bool test2(1);
   double thresh(1e-10);
   std::cout<<"Testing Jint_signflip_renormalize"<<std::endl;
   int l1,l2,l3,m1,m2,m3;

   //Test 1 Flip sign and keep value of m's
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      int cm1,cm2,cm3;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=(rand() % 10);
       m1=( -l1 + rand() % (2*l1+2));
       m2=( -l2 + rand() % (2*l3+2));
       m3=( -l3 + rand() % (2*l2+2));

       //Copy the values
       cm1=m1;
       cm2=m2;
       cm3=m3;

       //Send the parameters for sign flip
       Jint_signflip_renormalize(l1,l2,l3,&m1,&m2,&m3);

       //Check if the sign is flipped and the value conserved

       if( m1 == abs(cm1) && m2 == abs(cm2) && m3 == abs(cm3))
          test1*=1;
       else
          test1*=0;
   }
   //Test 2 Properly renormalized
   std::cout<<"2..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      double val;
      double check;
      int cm1,cm2,cm3;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       l3=(rand() % 5);
       m1=( -l1 + rand() % (2*l1+1));
       m2=( -l2 + rand() % (2*l2+1));
       m3=( -l3 + rand() % (2*l3+1));

       //Copy the values
       cm1=m1;
       cm2=m2;
       cm3=m3;

       //Send the parameters for sign flip
       val=Jint_signflip_renormalize(l1,l2,l3,&m1,&m2,&m3);

//       std::cout<<l1<<","<<l2<<","<<l3<<" - "<<cm1<<","<<cm2<<","<<cm3<<"+++++"<<val<<" ----- ";
       //Compute the renormalization constant using independent functions
       check=1;
       if(cm1 < 0)
           check*= pow( -1 , m1 ) * exp( lgamma(l1-m1+1) - lgamma(l1+m1+1) ); 
       if(cm2 < 0)
           check*= pow( -1 , m2 ) * exp( lgamma(l2-m2+1) - lgamma(l2+m2+1) );
       if(cm3 < 0)
           check*= pow( -1 , m3 ) * exp( lgamma(l3-m3+1) - lgamma(l3+m3+1) ); 

//       std::cout<<check<<std::endl;
       //Check if the renormalization constant is correct

       if( fabs(val - check) / (check + bool(check == 0)) <=thresh)
          test2*=1;
       else
          test2*=0;
   }
   if(test1 && test2)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      if(!test2)
         std::cout<<"Error 2...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_Jint_normalize()
{
   bool test1(1);
   double thresh(1e-10);
   std::cout<<"Testing Jint_normalize"<<std::endl;
   int l1,l2,l3,m1,m2,m3;
   //Test 1 Properly normalized
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       l3=(rand() % 5);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the normalization
//       std::cout<<"====="<<val<<std::endl;
       val=Jint_normalize(l1,l2,l3,m1,m2,m3);
//       std::cout<<"=====>>>"<<val<<std::endl;

//       std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
       //Compute the renormalization constant using independent functions
       check=1;
           check*= exp( - ( lgamma(l1-m1+1) - lgamma(l1+m1+1) ) / 2 ); 
           check*= exp( - ( lgamma(l2-m2+1) - lgamma(l2+m2+1) ) / 2 );
           check*= exp( - ( lgamma(l3-m3+1) - lgamma(l3+m3+1) ) / 2 ); 

//       std::cout<<check<<std::endl;
       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check == 0)) <=thresh)
          test1*=1;
       else
          test1*=0;
   }
   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_Jint_special_cases()
{
   bool test1(1);
   bool test2(1);
   bool test3(1);
   double thresh(1e-5);
   std::cout<<"Testing Jint_special_cases"<<std::endl;
   int l1,l2,l3,m1,m2,m3;


   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1+2*double(i)/double(nx);
   double dx(x[1]-x[0]);


   //Test 1 : Returns 0 if it is not a special case, or 1 if it is a special case
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      bool ret(NAN);
      double val(10);
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=(rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       ret=Jint_special_cases(l1,l2,l3,m1,m2,m3,&val); 
//       std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
       //Assess if the case is a special case 
       if((m3==m1+m2 || m3 == fabs(m1-m2) || m1>l1 || m2>l2 || m3>l3 || l1<0 || l2<0 || l3<0 || (l1+l2+l3+m1+m2+m3) != 0) && ret)
          test1*=1;
       else if(!(m3==m1+m2 || m3 == fabs(m1-m2) || m1>l1 || m2>l2 || m3>l3 || l1<0 || l2<0 || l3<0 || (l1+l2+l3+m1+m2+m3) != 0) && !ret)
          test1*=1;
       else
          test1*=1;
   }


   //Test 2 : The integral yields zero if the parity is odd
   std::cout<<"2..."<<std::endl;
   bool integrated(0);
   for(int i=0;i!=100;i++)
   {
      bool ret(NAN);
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       l3=(rand() % 5);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       if( ( l1 + l2 + l3 + m1 + m2 + m3 ) % 2 == 0 )
          l1++;

       //Send the parameters for computing the special cases
       ret=Jint_special_cases(l1,l2,l3,m1,m2,m3,&val); 

       if(ret)
       {
          integrated=1;
//       std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          //Compute the renormalization constant using independent functions
          check=0;
          for(int i=0;i!=nx;i++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i])*gsl_sf_legendre_Plm(l2,m2,x[i])*gsl_sf_legendre_Plm(l3,m3,x[i]);
//       std::cout<<check<<std::endl;

       //Check if the normalization constant is correct

       if( fabs(val - check) <=thresh)
          test2*=1;
       else
          test2*=0;
       }
   }
   if(!integrated)
   {
      test2*=0;
      std::cout<<" NO INTEGRATION DONE IN TEST 2 :"<<std::endl;
   }


   //Test 3 : The integral with m3=m1+m2
   std::cout<<"3..."<<std::endl;
   for(int i=0;i!=100;i++)
   {
      bool ret(NAN);
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       l3=(m1+m2 + rand() % 5);
       //m3=( rand() % (l3+1));
       m3=m1+m2;

       //Send the parameters for computing the special cases
       ret=Jint_special_cases(l1,l2,l3,m1,m2,m3,&val); 

       if(ret)
       {
          integrated=1;
//       std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          //Compute the renormalization constant using independent functions
          check=0;
          for(int i=0;i!=nx;i++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i])*gsl_sf_legendre_Plm(l2,m2,x[i])*gsl_sf_legendre_Plm(l3,m3,x[i]);
//       std::cout<<check<<std::endl;

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test3*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test3*=0;
       }
       }
   }
   if(!integrated)
   {
      test3*=0;
      std::cout<<" NO INTEGRATION DONE IN TEST 3 "<<std::endl;
   }

   if(test1 && test2 && test3)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      if(!test2)
         std::cout<<"Error 2...";
      if(!test3)
         std::cout<<"Error 3...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_ALP_integral()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing ALP_integral"<<std::endl;
   int l1,m1;

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1+2*double(i)/double(nx);
   double dx(x[1]-x[0]);

   //Test 1 : The ALP integral is correct
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       m1=( rand() % (l1+1));

       //Send the parameters for computing the special cases
       val=ALP_integral(l1,m1); 

     //std::cout<<l1<<" - "<<m1<<"+++++"<<val<<" ----- ";
       //Compute the renormalization constant using independent functions
       check=0;
       for(int i=0;i!=nx;i++)
          check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i]);
    // std::cout<<check<<std::endl;

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<" - "<<m1<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_two_ALP_integral()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing two_ALP_integral"<<std::endl;
   int l1,l2,m1,m2;


   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1+2*double(i)/double(nx);
   double dx(x[1]-x[0]);

   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));

       //Send the parameters for computing the special cases
       val=two_ALP_integral(l1,l2,m1,m2); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int i=0;i!=nx;i++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i])*gsl_sf_legendre_Plm(l2,m2,x[i]);

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<" - "<<m1<<","<<m2<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_three_ALP_J_integral()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing three_ALP_J_integral"<<std::endl;
   int l1,l2,l3,m1,m2,m3;


   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1+2*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=( rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=three_ALP_J_integral(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int i=0;i!=nx;i++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i])*gsl_sf_legendre_Plm(l2,m2,x[i])*gsl_sf_legendre_Plm(l3,m3,x[i]);

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}

//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_I_m1_integral()
{
   std::cout<<"Testing I_m1_integral"<<std::endl;
   int m1,m2,m3;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);
   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral give good values
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);
       m3=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*sin(x[i])*phi(m1,x[i])*phi(m2,x[i])*phi(m3,x[i]);
       
       //generate function result
       val=I_m1_integral(m1,m2,m3);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<","<<m3<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_I_p1_integral()
{
   std::cout<<"Testing I_p1_integral"<<std::endl;
   int m1,m2,m3;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);
   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral give good values
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);
       m3=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*cos(x[i])*phi(m1,x[i])*phi(m2,x[i])*phi(m3,x[i]);
       
       //generate function result
       val=I_p1_integral(m1,m2,m3);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<","<<m3<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_I_m1_D_integral()
{
   std::cout<<"Testing I_m1_D_integral"<<std::endl;
   int m1,m2,m3;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);
   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral give good values
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=1000;j++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);
       m3=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*sin(x[i])*phi(m1,x[i])*(-m2*phi(-m2,x[i])*phi(m3,x[i])-m3*phi(m2,x[i])*phi(-m3,x[i]));
       
       //generate function result
       val=I_m1_D_integral(m1,m2,m3);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<","<<m3<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_I_p1_D_integral()
{
   std::cout<<"Testing I_p1_D_integral"<<std::endl;
   int m1,m2,m3;
   double sum;
   double val;
   double thresh(1e-6);

   bool test1(1);
   //set up vector tor epresent function
   int nx=100000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=2*acos(-1)*i/(nx);
   double dx(x[1]-x[0]);


   //Initialize random variable generation
   srand (time(0));

   //test 1 : the integral give good values
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=1000;j++)
   {
       m1=-10 + (rand() % 21);
       m2=-10 + (rand() % 21);
       m3=-10 + (rand() % 21);

       //Generate numerical value
       sum=0;
       for(int i=0;i!=nx;i++)
          sum+=dx*cos(x[i])*phi(m1,x[i])*(-m2*phi(-m2,x[i])*phi(m3,x[i])-m3*phi(m2,x[i])*phi(-m3,x[i]));
       
       //generate function result
       val=I_p1_D_integral(m1,m2,m3);

       //compare numerical value to function result
       if(val==0 && fabs(sum)<=thresh)
          test1*=1;
       else if(val!=0 && fabs(val-sum)<=thresh)
          test1*=1;
       else
       {
          test1*=0;
          std::cout<<m1<<","<<m2<<","<<m3<<" : "<<val<<" ; "<<sum<<std::endl;
       }
   }

   if(test1 )
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_J_int_m2()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing J_int_m2"<<std::endl;
   int l1,l2,l3,m1,m2,m3;

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-.99999+1.99999*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=( rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=J_int_m2(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int i=0;i!=nx;i++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[i])*gsl_sf_legendre_Plm(l2,m2,x[i])*gsl_sf_legendre_Plm(l3,m3,x[i])/sqrt(1-x[i]*x[i]);

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}

//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_J_int_m1()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing J_int_m1"<<std::endl;
   int l1,l2,l3,m1,m2,m3;

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1.+2.*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=100;i++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=( rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=J_int_m1(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int xx=0;xx!=nx;xx++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[xx])*gsl_sf_legendre_Plm(l2,m2,x[xx])*gsl_sf_legendre_Plm(l3,m3,x[xx])*sqrt(1-x[xx]*x[xx]);

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_J_int_p1()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing J_int_p1"<<std::endl;
   int l1,l2,l3,m1,m2,m3;

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1.+2.*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=100;i++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 10);
       l2=(rand() % 10);
       l3=( rand() % 10);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=J_int_p1(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int xx=0;xx!=nx;xx++)
             check+=dx*gsl_sf_legendre_Plm(l1,m1,x[xx])*gsl_sf_legendre_Plm(l2,m2,x[xx])*gsl_sf_legendre_Plm(l3,m3,x[xx])*x[xx];

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_J_int_m1_D()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing J_int_m1_D"<<std::endl;
   int l1,l2,l3,m1,m2,m3;
   double P1,P2,P3;
   double DP2,DP3;

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-.999999+1.999999*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=100;i++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       l3=( rand() % 5);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=J_int_m1_D(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int xx=0;xx!=nx;xx++)
          {
             P1=(gsl_sf_legendre_Plm(l1,m1,x[xx]));
             P2=(gsl_sf_legendre_Plm(l2,m2,x[xx]));
             P3=(gsl_sf_legendre_Plm(l3,m3,x[xx]));
             if(l2==0)
                DP2=0;
             else if(l2>0 && m2==0)
                DP2=gsl_sf_legendre_Plm(l2,1,x[xx]);
             else if(l2==m2)
                DP2=-0.5*((l2+m2)*(l2-m2+1)*gsl_sf_legendre_Plm(l2,m2-1,x[xx]));
             else
                DP2=-0.5*((l2+m2)*(l2-m2+1)*gsl_sf_legendre_Plm(l2,m2-1,x[xx])-gsl_sf_legendre_Plm(l2,m2+1,x[xx]));
             if(l3==0)
                DP3=0;
             else if(l3>0 && m3==0)
                DP3=gsl_sf_legendre_Plm(l3,1,x[xx]);
             else if(l3==m3)
                DP3=0.5*((l3+m3)*(l3-m3+1)*gsl_sf_legendre_Plm(l3,m3-1,x[xx]));
             else
                DP3=-0.5*((l3+m3)*(l3-m3+1)*gsl_sf_legendre_Plm(l3,m3-1,x[xx])-gsl_sf_legendre_Plm(l3,m3+1,x[xx]));

             check+=dx*(P1*(DP2*P3+DP3*P2));
          }

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_J_int_p1_D()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing J_int_p1_D"<<std::endl;
   int l1,l2,l3,m1,m2,m3;
   double P1,P2,P3;
   double DP2,DP3;

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=5000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=-1.+2.*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=150;i++)
   {
      double val(10);
      double check;
      //generate a random set of values
       l1=(rand() % 5);
       l2=(rand() % 5);
       l3=( rand() % 5);
       m1=( rand() % (l1+1));
       m2=( rand() % (l2+1));
       m3=( rand() % (l3+1));

       //Send the parameters for computing the special cases
       val=J_int_p1_D(l1,l2,l3,m1,m2,m3); 

          //Compute the renormalization constant using independent functions
          check=0;
          for(int xx=0;xx!=nx;xx++)
          {
             P1=(gsl_sf_legendre_Plm(l1,m1,x[xx]));
             P2=(gsl_sf_legendre_Plm(l2,m2,x[xx]));
             P3=(gsl_sf_legendre_Plm(l3,m3,x[xx]));
             if(l2==0)
                DP2=0;
             else if(l2>0 && m2==0)
                DP2=gsl_sf_legendre_Plm(l2,1,x[xx]);
             else if(l2==m2)
                DP2=-0.5*((l2+m2)*(l2-m2+1)*gsl_sf_legendre_Plm(l2,m2-1,x[xx]));
             else
                DP2=-0.5*((l2+m2)*(l2-m2+1)*gsl_sf_legendre_Plm(l2,m2-1,x[xx])-gsl_sf_legendre_Plm(l2,m2+1,x[xx]));
             if(l3==0)
                DP3=0;
             else if(l3>0 && m3==0)
                DP3=gsl_sf_legendre_Plm(l3,1,x[xx]);
             else if(l3==m3)
                DP3=0.5*((l3+m3)*(l3-m3+1)*gsl_sf_legendre_Plm(l3,m3-1,x[xx]));
             else
                DP3=-0.5*((l3+m3)*(l3-m3+1)*gsl_sf_legendre_Plm(l3,m3-1,x[xx])-gsl_sf_legendre_Plm(l3,m3+1,x[xx]));

             check+=x[xx]*dx*(P1*(DP2*P3+DP3*P2));
          }

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
          test1*=1;
       else
       {
          std::cout<<l1<<","<<l2<<","<<l3<<" - "<<m1<<","<<m2<<","<<m3<<"+++++"<<val<<" ----- ";
          std::cout<<check<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_B_coeff()
{
   bool test1(1);
   double thresh(1e-5);
   std::cout<<"Testing B_coeff"<<std::endl;
   int l1,m1,m2;
   int m;
   double thet,phi;
   std::vector<double> B_val;
   std::complex<double> val;

   //Initialize random variable generation
   srand (time(0));

   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int j=0;j!=100;j++)
   {
      std::complex<double> check;
      //generate a random set of values
       l1=(rand() % 10);
       m1=( rand() % ( 2 * l1 + 1 ) ) - l1;
       m2=( rand() % 3 ) - 1;
       thet=( rand() % (1000) * acos(-1) / 1000 );
       phi=( rand() % (1000) * 2 * acos(-1) / 1000 );

       //Send the parameters for computing the special cases
       B_coeff(l1,m1,m2,&B_val); 

      //Check for a definite relation using independent functions
      check=pow(-1,m1*bool(m1<0))*exp(std::complex<double>(0,1)*m1*phi)*gsl_sf_legendre_sphPlm(l1,abs(m1),cos(thet))
         *pow(-1,m2*bool(m2<0))*exp(std::complex<double>(0,1)*m2*phi)*gsl_sf_legendre_sphPlm(1,abs(m2),cos(thet));

      m=m1+m2;
      if((l1-1)>=abs(m) && l1 > 0)
      {
         val=B_val.at(0)*pow(-1,m*bool(m<0))*exp(std::complex<double>(0,1)*m*phi)*gsl_sf_legendre_sphPlm(l1-1,abs(m),cos(thet))
            +B_val.at(1)*pow(-1,m*bool(m<0))*exp(std::complex<double>(0,1)*m*phi)*gsl_sf_legendre_sphPlm(l1+1,abs(m),cos(thet));
      }

      else
         val=B_val.at(1)*pow(-1,m*bool(m<0))*exp(std::complex<double>(0,1)*m*phi)*gsl_sf_legendre_sphPlm(l1+1,abs(m),cos(thet));

       //Check if the relation holds

       if( abs(check-val)<=thresh)
       {
          test1*=1;
       }
       else
       {
          std::cout<<l1<<" - "<<m1<<","<<m2<<" ----- ";
          std::cout<<check<<"====="<<val<<std::endl;
          test1*=0;
       }
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
