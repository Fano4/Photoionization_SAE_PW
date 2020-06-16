#include "test_utils.cpp"
#include "test_angular_int_aux.cpp"
#include "test_files_utils.cpp"
#include "test_molec_integ.cpp"
#include "test_poly_ovlp_one_elec_ope.cpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <omp.h>

void two_file_overlap();

int main(int argc,char* argv[])
{
   omp_set_num_threads(8);
   if(argc<2)
   {
      // report version
      std::cout<<argv[0]<<" : Testing PhotoICE code."<<std::endl;

   }

   two_file_overlap();
   exit(EXIT_SUCCESS);
//   test_fact_prime_decomposer();
//   two_azim_integ_test();
//   three_azim_integ_test();
//   test_I_m1_integral();
//   test_I_p1_integral();
//   test_I_m1_D_integral();
//   test_Jint_sort_indices();
//   test_Jint_signflip_renormalize();
//   test_Jint_normalize();
//   test_Jint_special_cases();
//   test_ALP_integral();
//   test_two_ALP_integral();
//   test_three_ALP_J_integral();
//   test_J_int_m2();
//   test_J_int_m1();
//   test_J_int_p1();
//   test_J_int_m1_D();
//   test_J_int_p1_D();
//   test_search();
   test_prim_ovlp();
   test_molp_sym_parser();
   test_molp_method_parser();
   test_molp_wf_parser();
   test_molp_cas_reader();
   test_molp_basis_parser();
   test_molp_lcao_parser();
   test_molp_ci_parser();
   test_prim_radial_ovlp();
   test_ao_ovlp();
   test_mo_ovlp();
   test_slater_ovlp();
   test_es_ovlp(); 
   return 0;
}
void two_file_overlap()
{

   std::string out_file("/home/users/stephan/cayo_inputs/neut_es_ovlp.txt");
   std::string coord_file("/home/users/stephan/cayo_inputs/ch4_path_file.out");
   std::string fileroot("");
   std::stringstream sstr;

   std::ifstream coord_istr;
   std::vector<std::string> x;

   int n_states(1);
   std::vector<double> ES_MO;

   coord_istr.open(coord_file.c_str());

   if(!coord_istr.is_open())
   {
      std::cout<<"could not open coord file"<<std::endl;
   }
   std::string tmp_str;
   while(!coord_istr.eof())
   {
      coord_istr>>tmp_str;
      x.push_back(tmp_str);
   }
   coord_istr.close();
   int N_geom(x.size()-1);
  
   std::cout<<"Starting overlap between different geometries. There are "<<N_geom<<" geometries to compare"<<std::endl;
   std::cout<<"We have to evaluate "<<N_geom-1<<"overlaps"<<std::endl;

   std::ofstream output;
   output.open(out_file.c_str(),std::ios_base::trunc);
   //add the first point to obtain the original grid dimensionality
   output<<std::fixed<<std::setprecision(8);//<<x[0];
   for(int i=0;i!=n_states;i++)
      output<<std::setw(12)<<1<<" ";
   output<<std::endl;
   for(int i=0;i!=N_geom-1;i++)
   {
      std::cout<<"Now reading file "<<x[i]<<" and "<<x[i+1]<<std::endl;
      sstr.str(x[i].c_str());
//      sstr<<fileroot<<x[i]<<".out";
      std::string file1_loc(sstr.str());
      sstr.str(x[i+1].c_str());
//      sstr<<fileroot<<x[i+1]<<".out";
      std::string file2_loc(sstr.str());


      test_es_ovlp_twogeoms(file1_loc,file2_loc,&ES_MO);

      output<<std::fixed<<std::setprecision(8);//<<x[i+1];
      for(int j=0;j!=n_states;j++)
         output<<std::setw(12)<<ES_MO.at(j*n_states+j)<<" ";
      output<<std::endl;
   }

   output.close();
}
