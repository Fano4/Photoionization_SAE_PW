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
bool test_es_ovlp()
{
   using namespace std;

   bool test1(0);

   cout<<"testing es_ovlp...";

   string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   //set up the environment to compute the overlap. Computing the slater overlap requires the following data

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
   vector<int> n_occ;
   vector<int> n_closed;
   vector<int> n_frozen;
   vector<int> csf_mo;
   vector<int> csf_spin;
   vector<double> ci_coeff_sym;
   vector<int> ci_num_sym;
   vector<int> n_elec;
   vector<int> sym;
   vector<int> spin;
   vector<int> charge;
   vector<int> n_states;


   //Get the basic data from the input file
   n_sym=molp_sym_parser(input_file);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn,&yn,&zn,input_file);
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l,&ml,&cont_zeta,&cont_coeff,input_file);
   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file);
   molp_wf_parser(method_index,&n_elec,&sym,&spin,&charge,&n_states,input_file);
   molp_lcao_parser(method_index,&lcao_coeff_sym,input_file);
   molp_ci_parser(method_index,&csf_mo,&csf_spin,&ci_coeff_sym,&ci_num_sym,input_file);


   //Transform geometry data to a single vector
   for(int i=0;i<xn.size();i++)
   {
      cart_r.push_back(xn.at(i));
      cart_r.push_back(yn.at(i));
      cart_r.push_back(zn.at(i));
   }



   //Initiate the lcao_coeff array in single block
   int basis_size_tot;
   int n_occ_tot;
   int n_states_tot;
   int ci_num_tot;
   vector<double> lcao_coeff;
   vector<double> ci_coeff;

   //Transform lcao coeff matrix and ci vector in single block
   sym_to_nosym_mat_trans(n_sym,n_occ,basis_size,lcao_coeff_sym,&n_occ_tot,&basis_size_tot,&lcao_coeff);
   sym_to_nosym_mat_trans(n_sym,ci_num_sym,n_states,ci_coeff_sym,&ci_num_tot,&n_states_tot,&ci_coeff);

   //Compute ao overlap matrix
   vector<double> S;
   ao_ovlp(cart_r,cart_r,nuc_bas_func,nuc_bas_func,cont_num,cont_num,cont_zeta,cont_zeta,cont_coeff,cont_coeff,l,l,ml,ml,&S);
   std::cout<<"The S matrix contains "<<S.size()<<" elements."<<std::endl;
   
   //Compute mo_overlap matrix
   vector<double> MO_S;
   MO_ovlp(S,lcao_coeff,lcao_coeff,&MO_S);
   std::cout<<"The MO S matrix contains "<<MO_S.size()<<" elements."<<std::endl;

   //Compute the overlap between the Slater determinants
   vector<double> CSF_S;
   slater_ovlp(n_elec.at(0),ci_num_tot,ci_num_tot,csf_mo,csf_mo,csf_spin,csf_spin,n_occ_tot,n_occ_tot,MO_S,&CSF_S);

   //Compute the overlap between electronic states
   vector<double> ES_S;
   ES_ovlp(CSF_S,ci_num_tot,ci_num_tot,ci_coeff,ci_coeff,n_states_tot,n_states_tot,&ES_S);

   for(int i=0;i!=n_states_tot;i++)
   {
      for(int j=0;j!=n_states_tot;j++)
      {
         std::cout<<i<<" , "<<j<<" - "<<ES_S.at(i*n_states_tot+j)<<std::endl;
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
bool test_es_ovlp_twogeoms(std::string input_file_a,std::string input_file_b,std::vector<double>* ES_S)
{
   using namespace std;

   bool test1(0);

   cout<<"testing es_ovlp in multiple files...";

   //set up the environment to compute the overlap. Computing the slater overlap requires the following data

   int method_index(0);
   int num_of_nucl;
   int n_sym;
   vector<double> cart_r_a;
   vector<double> xn_a;
   vector<double> yn_a;
   vector<double> zn_a;
   vector<double> cart_r_b;
   vector<double> xn_b;
   vector<double> yn_b;
   vector<double> zn_b;
   vector<double> cont_zeta;
   vector<double> cont_coeff;
   vector<double> lcao_coeff_sym_a;
   vector<double> lcao_coeff_sym_b;
   vector<int> basis_size;
   vector<int> cont_num;
   vector<int> cont_val;
   vector<int> nuc_bas_func;
   vector<int> Z_nucl;
   vector<unsigned int> l;
   vector<int> ml;
   vector<int> n_occ;
   vector<int> n_closed;
   vector<int> n_frozen;
   vector<int> csf_mo_a;
   vector<int> csf_mo_b;
   vector<int> csf_spin_a;
   vector<int> csf_spin_b;
   vector<double> ci_coeff_sym_a;
   vector<double> ci_coeff_sym_b;
   vector<int> ci_num_sym_a;
   vector<int> ci_num_sym_b;
   vector<int> n_elec;
   vector<int> sym;
   vector<int> spin;
   vector<int> charge;
   vector<int> n_states;


   //Get the basic data from the input file
   n_sym=molp_sym_parser(input_file_a);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn_a,&yn_a,&zn_a,input_file_a);
   molp_geom_parser(&num_of_nucl,&Z_nucl,&xn_b,&yn_b,&zn_b,input_file_b);
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l,&ml,&cont_zeta,&cont_coeff,input_file_a);
   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file_a);
   molp_wf_parser(method_index,&n_elec,&sym,&spin,&charge,&n_states,input_file_a);
   molp_lcao_parser(method_index,&lcao_coeff_sym_a,input_file_a);
   molp_lcao_parser(method_index,&lcao_coeff_sym_b,input_file_b);
   molp_ci_parser(method_index,&csf_mo_a,&csf_spin_a,&ci_coeff_sym_a,&ci_num_sym_a,input_file_a);
   molp_ci_parser(method_index,&csf_mo_b,&csf_spin_b,&ci_coeff_sym_b,&ci_num_sym_b,input_file_b);


   //Transform geometry data to a single vector
   std::cout<<xn_a.size()<<","<<xn_b.size()<<std::endl;
   for(int i=0;i<xn_a.size();i++)
   {
      cart_r_a.push_back(xn_a.at(i));
      cart_r_a.push_back(yn_a.at(i));
      cart_r_a.push_back(zn_a.at(i));
      cart_r_b.push_back(xn_b.at(i));
      cart_r_b.push_back(yn_b.at(i));
      cart_r_b.push_back(zn_b.at(i));
   }



   //Initiate the lcao_coeff array in single block
   int basis_size_tot;
   int n_occ_tot;
   int n_states_tot;
   int ci_num_tot_a;
   int ci_num_tot_b;
   vector<double> lcao_coeff_a;
   vector<double> lcao_coeff_b;
   vector<double> ci_coeff_a;
   vector<double> ci_coeff_b;

   //Transform lcao coeff matrix and ci vector in single block
   sym_to_nosym_mat_trans(n_sym,n_occ,basis_size,lcao_coeff_sym_a,&n_occ_tot,&basis_size_tot,&lcao_coeff_a);
   sym_to_nosym_mat_trans(n_sym,ci_num_sym_a,n_states,ci_coeff_sym_a,&ci_num_tot_a,&n_states_tot,&ci_coeff_a);
   sym_to_nosym_mat_trans(n_sym,n_occ,basis_size,lcao_coeff_sym_b,&n_occ_tot,&basis_size_tot,&lcao_coeff_b);
   sym_to_nosym_mat_trans(n_sym,ci_num_sym_b,n_states,ci_coeff_sym_b,&ci_num_tot_b,&n_states_tot,&ci_coeff_b);

   //Compute ao overlap matrix
   vector<double> S;
   ao_ovlp(cart_r_a,cart_r_b,nuc_bas_func,nuc_bas_func,cont_num,cont_num,cont_zeta,cont_zeta,cont_coeff,cont_coeff,l,l,ml,ml,&S);
   std::cout<<"The S matrix contains "<<S.size()<<" elements."<<std::endl;
   
   //Compute mo_overlap matrix
   vector<double> MO_S;
   MO_ovlp(S,lcao_coeff_a,lcao_coeff_b,&MO_S);
   std::cout<<"The MO S matrix contains "<<MO_S.size()<<" elements."<<std::endl;

   //Compute the overlap between the Slater determinants
   vector<double> CSF_S;
   slater_ovlp(n_elec.at(0),ci_num_tot_a,ci_num_tot_b,csf_mo_a,csf_mo_b,csf_spin_a,csf_spin_b,n_occ_tot,n_occ_tot,MO_S,&CSF_S);

   //Compute the overlap between electronic states
   ES_ovlp(CSF_S,ci_num_tot_a,ci_num_tot_b,ci_coeff_a,ci_coeff_b,n_states_tot,n_states_tot,ES_S);


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
