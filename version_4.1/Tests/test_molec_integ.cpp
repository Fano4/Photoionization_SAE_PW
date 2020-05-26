#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
//#include "filesutils.h"
#include "utilities.h"
#include "mathfunctions.h"
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

   //set up vector to represent function
   int nx=1000000;
   double*x=new double[nx];
   for(int i=0;i!=nx;i++)
      x[i]=20*double(i)/double(nx);
   double dx(x[1]-x[0]);
   //Test 1 : The integral yields good results
   std::cout<<"1..."<<std::endl;
   for(int i=0;i!=25;i++)
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
//          std::cout<<la<<","<<lb<<","<<l<<","<<zeta_a<<","<<zeta_b<<","<<r<<"+++++"<<val<<" ----- ";
//          std::cout<<check<<"   PASSED"<<std::endl;
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
bool test_prim_ovlp()
{

   using namespace std; 

   bool test1(0);

   cout<<"testing prim_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data




   vector<double> ra;
   vector<double> rb;

   ra.push_back(0);
   ra.push_back(0);
   ra.push_back(-1.);
   rb.push_back(0);
   rb.push_back(0);
   rb.push_back(2.);
   for(int i=0;i!=10000;i++)
      std::cout<<prim_ovlp(ra,rb,0.3,1000-0.1*i,0,0,0,0)<<std::endl;


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

}

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
   vector<unsigned int> l;
   vector<int> ml;

   n_sym=molp_sym_parser(input_file);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn,&yn,&zn,input_file);
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l,&ml,&cont_zeta,&cont_coeff,input_file);

   for(int i=0;i<Z_nucl.size();i++)
   {
      cart_r.push_back(xn.at(i));
      cart_r.push_back(yn.at(i));
      cart_r.push_back(zn.at(i));
   }

   int basis_size_tot(0);
   for(int i=0;i!=n_sym;i++)
      basis_size_tot+=basis_size[i];

   vector<double> S;
   ao_ovlp(cart_r,cart_r,nuc_bas_func,nuc_bas_func,cont_num,cont_num,cont_zeta,cont_zeta,cont_coeff,cont_coeff,l,l,ml,ml,&S);


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
         countb=memb;
         for(int j=0;j!=basis_size[mm];j++)
         {
            if(j%10 == 0 )
               cout<<endl;
            cout<<std::scientific<<std::setw(14)<<S.at(counta*basis_size_tot+countb);
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
bool test_mo_ovlp()
{
   using namespace std; 

   bool test1(0);

   cout<<"testing mo_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the ao_basis overlap requires the geometry and the basis data

   int method_index(0);
   int num_of_nucl;
   int n_sym;
   vector<double> cart_r;
   vector<double> xn;
   vector<double> yn;
   vector<double> zn;
   vector<double> cont_zeta;
   vector<double> cont_coeff;
   vector<double> lcao_coeff_sym;
   vector<int> basis_size;
   vector<int> cont_num;
   vector<int> cont_val;
   vector<int> nuc_bas_func;
   vector<int> Z_nucl;
   vector<unsigned int> l;
   vector<int> ml;
   std::vector<int> n_occ;
   std::vector<int> n_closed;
   std::vector<int> n_frozen;

   n_sym=molp_sym_parser(input_file);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn,&yn,&zn,input_file);
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l,&ml,&cont_zeta,&cont_coeff,input_file);
   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file);

   for(int i=0;i<xn.size();i++)
   {
      cart_r.push_back(xn.at(i));
      cart_r.push_back(yn.at(i));
      cart_r.push_back(zn.at(i));
   }

   molp_lcao_parser(method_index,&lcao_coeff_sym,input_file);

   int basis_size_tot(0);
   int n_occ_tot(0);

   vector<double> lcao_coeff(n_occ_tot*basis_size_tot);
   sym_to_nosym_mat_trans(n_sym,n_occ,basis_size,lcao_coeff_sym,&n_occ_tot,&basis_size_tot,&lcao_coeff);

   vector<double> S;
   ao_ovlp(cart_r,cart_r,nuc_bas_func,nuc_bas_func,cont_num,cont_num,cont_zeta,cont_zeta,cont_coeff,cont_coeff,l,l,ml,ml,&S);
   std::cout<<"The S matrix contains "<<S.size()<<" elements."<<std::endl;
   vector<double> MO_S;
   MO_ovlp(S,lcao_coeff,lcao_coeff,&MO_S);
   std::cout<<"The MO S matrix contains "<<MO_S.size()<<" elements."<<std::endl;

   int counta(0);
   int countb(0);
   int memb;

   std::cout<<"1"<<std::endl;
   for(int mm=0;mm<n_sym;mm++)
   {
      std::cout<<std::endl<<"SYMMETRY BLOCK "<<mm+1<<std::endl;
      memb=countb;

      for(int i=0;i!=n_occ[mm];i++)
      {
         countb=memb;
         for(int j=0;j!=n_occ[mm];j++)
         {
            if(j%10 == 0 )
               cout<<endl;
            cout<<std::scientific<<std::setw(14)<<MO_S.at(counta*n_occ_tot+countb);
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
