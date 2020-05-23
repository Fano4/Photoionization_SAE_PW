#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
//#include "filesutils.h"
#include "mathfunctions.h"
/*
bool test_bessel_gaussian_poly_integral()
{

   using namespace std; 

   bool test1(1);
   unsigned int l1,l2;
   double zeta;
   double r;
   double thresh(1e-10);

   cout<<"testing bessel_gaussian_poly_integral...";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data

   //Initialize random variable generation
   srand (time(0));

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=10*double(i)/double(nx);
   double dx(x[1]-x[0]);
   
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
   {
      double val(10);
      double check;

      l1 = (rand() % 6);
      l2 = (rand() % (l1+1));
      l1 = 2*l1-l2;
      zeta = 0.1 + 100 * double( rand() % 1000 ) / 1000;
      r = 0.5 + 3. * double( rand() % 1000 ) / 1000;
      //Send the parameters for computing 
      val=bessel_gaussian_poly_integral(l1,l2,zeta,r);
      //Compute using independent functions
      check=0;
      for(int xx=0;xx!=nx;xx++)
         check+=dx*pow(x[xx],2+l1)*exp(-zeta*x[xx]*x[xx])*gsl_sf_bessel_jl(l2,x[xx]*r);

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
       {
          //std::cout<<l1<<","<<l2<<","<<zeta<<","<<r<<"+++++"<<val<<" ----- ";
          //std::cout<<check<<"  OK "<<std::endl;
          test1*=1;
       }
       else
       {
          std::cout<<l1<<","<<l2<<","<<zeta<<","<<r<<"+++++"<<val<<" ----- ";
          std::cout<<check<<"  FAILED "<<std::endl;
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
      {
         std::cout<<"Error 1...";
      }
      std::cout<<std::endl;
      return 0;
   }

}*/
bool test_prim_radial_ovlp()
{

   using namespace std; 

   bool test1(0);
   double thresh(1e-8);

   cout<<"testing prim_radial_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data


   double zeta_a;
   double zeta_b;
   unsigned int la;
   unsigned int lb;
   unsigned int l;
   double r;

   //set up vector tor epresent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=10*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=100;i++)
   {
      double val(10);
      double a,b;
      double check;

      la = (rand() % 6);
      lb = (rand() % 6);
      l = abs(int(la-lb))+(rand() % (la+lb+1));
      if((la+lb+l)%2 != 0)
      {
         i--;
         continue;
      }
      zeta_a = 0.001 + 10. * double ( rand() % 1000 ) / 1000.;
      zeta_b = 0.001 + 10. * double ( rand() % 1000 ) / 1000.;
      r = 0.5 + 3. * double( rand() % 1000 ) / 1000.;
      //Send the parameters for computing 
      val=prim_radial_ovlp(la,lb,l,zeta_a,zeta_b,r);
      //Compute using independent functions
      check=0;
      a=2+la+lb;
      b=(zeta_a+zeta_b)/(4*zeta_b*zeta_a);
      for(int xx=0;xx!=nx;xx++)
         check+=dx*pow(x[xx],a)*exp(-b*x[xx]*x[xx])*gsl_sf_bessel_jl(l,x[xx]*r)/ ( pow(2.*zeta_a,1.5+la) * pow(2.*zeta_b, 1.5 + lb) ) ;

       //Check if the normalization constant is correct

       if( fabs(val - check) / (check + bool(check <= thresh)) <=thresh)
       {
          std::cout<<la<<","<<lb<<","<<l<<","<<zeta_a<<","<<zeta_b<<","<<r<<"+++++"<<val<<" ----- ";
          std::cout<<check<<"   PASSED"<<std::endl;
          test1*=1;
       }
       else
       {
          std::cout<<la<<","<<lb<<","<<l<<","<<zeta_a<<","<<zeta_b<<","<<r<<"+++++"<<val<<" ----- ";
          std::cout<<check<<"   FAILED : relative error of "<<fabs(val - check) / (check + bool(check <= thresh))<<std::endl;
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
      {
         std::cout<<"Error 1...";
      }
      std::cout<<std::endl;
      return 0;
   }

}
/*
bool test_prim_ovlp()
{

   using namespace std; 

   bool test1(0);

   cout<<"testing prim_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data




   vector<double> ra;
   vector<double> rb;



   for(int i=0;i!=basis_size_tot;i++)
   {
      ra.clear();
      ra.push_back(0);
      ra.push_back(0);
      ra.push_back(-1);
      for(int j=0;j!=basis_size_tot;j++)
      {
         rb.clear();
         rb.push_back(0);
         rb.push_back(0);
         rb.push_back(1);

         cout<<setw(12)<<prim_ovlp(ra.data(),rb.data(),zeta_a,zeta_b,0,0,0,0);
      }
   }


   test1=1;

   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
      {
         std::cout<<"Error 1...";
      }
      std::cout<<std::endl;
      return 0;
   }

}*/

bool test_ao_ovlp()
{

   using namespace std; 

   bool test1(0);

   cout<<"testing ao_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data

   int num_of_nucl;
   int n_sym;
   vector<double> cart_r;
   vector<double> xn;
   vector<double> yn;
   vector<double> zn;
   vector<double> cont_zeta;
   vector<double> cont_coeff;
   vector<int> basis_size;
   vector<int> cont_num;
   vector<int> cont_val;
   vector<int> nuc_bas_func;
   vector<int> Z_nucl;
   vector<int> l;
   vector<int> ml;

   n_sym=molp_sym_parser(input_file);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn,&yn,&zn,input_file);
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l,&ml,&cont_zeta,&cont_coeff,input_file);

   int basis_size_tot(0);
   for(int i=0;i!=n_sym;i++)
      basis_size_tot+=basis_size[i];

   vector<int> ao_prim_index;
   vector<double>* cont_zeta_vec=new vector<double> [basis_size_tot];
   vector<double>* cont_coeff_vec=new vector<double> [basis_size_tot];

   int count(0);
   for(int i=0;i!=basis_size_tot;i++)
   {
      ao_prim_index.push_back(count);
      for(int j=0;j!=cont_num[i];j++)
      {
         cont_zeta_vec[i].push_back(cont_zeta.at(count));
         cont_coeff_vec[i].push_back(cont_coeff.at(count));
         count++;
      }
   }

   vector<double> ra;
   vector<double> rb;
   int counta(0);
   int countb(0);
   int memb;

   std::cout<<"1"<<std::endl;
   for(int mm=0;mm<n_sym;mm++)
   {
      std::cout<<std::endl<<"SYMMETRY BLOCK "<<mm+1<<std::endl;
      memb=countb;

      for(int i=0;i!=basis_size[mm];i++)
      {
         ra.clear();
         ra.push_back(xn.at(nuc_bas_func.at(counta)));
         ra.push_back(yn.at(nuc_bas_func.at(counta)));
         ra.push_back(zn.at(nuc_bas_func.at(counta)));

         
         countb=memb;
         for(int j=0;j!=basis_size[mm];j++)
         {
            rb.clear();
            rb.push_back(xn.at(nuc_bas_func.at(countb)));
            rb.push_back(yn.at(nuc_bas_func.at(countb)));
            rb.push_back(zn.at(nuc_bas_func.at(countb)));

            if(j%10 == 0 )
               cout<<endl;

//         ao_ovlp(ra,rb,cont_zeta_vec[i],cont_zeta_vec[j],cont_coeff_vec[i],cont_coeff_vec[j],l[i],l[j],ml[i],ml[j]);
            cout<<std::scientific<<std::setw(14)<<ao_ovlp(ra,rb,cont_zeta_vec[counta],cont_zeta_vec[countb],cont_coeff_vec[counta],cont_coeff_vec[countb],l[counta],l[countb],ml[counta],ml[countb]);
            countb++;
         }cout<<endl;
         counta++;
      }

   }


   test1=1;

   if(test1)
   {
      std::cout<<std::defaultfloat<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
      {
         std::cout<<"Error 1...";
      }
      std::cout<<std::endl;
      return 0;
   }

}
