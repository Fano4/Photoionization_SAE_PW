#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include "mathfunctions.h"

bool test_slater_ovlp()
{
   using namespace std;

   bool test1(0);

   cout<<"testing slater_ovlp...";

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

   for(int i=0;i!=ci_num_tot;i++)
   {
      vector<int> mo_a;
      vector<int> spin_a;
      for(int n=0;n!=n_elec.at(0);n++)
      {
         mo_a.push_back(csf_mo.at(i*n_elec.at(0)+n));
         spin_a.push_back(csf_spin.at(i*n_elec.at(0)+n));
      }
      for(int j=0;j!=ci_num_tot;j++)
      {
         vector<int> mo_b;
         vector<int> spin_b;
         for(int n=0;n!=n_elec.at(0);n++)
         {
            mo_b.push_back(csf_mo.at(j*n_elec.at(0)+n));
            spin_b.push_back(csf_spin.at(j*n_elec.at(0)+n));
         }
         std::cout<<i<<" - "<<j<<" : "<<slater_ovlp(mo_a,mo_b,spin_a,spin_b,MO_S)<<std::endl;
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
