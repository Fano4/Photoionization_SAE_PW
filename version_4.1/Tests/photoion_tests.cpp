#include "test_utils.cpp"
#include "test_angular_int_aux.cpp"
#include "test_files_utils.cpp"
#include "test_molec_integ.cpp"
#include <iostream>

int main(int argc,char* argv[])
{
   if(argc<2)
   {
      // report version
      std::cout<<argv[0]<<" : Testing PhotoICE code."<<std::endl;

   }
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
   test_molp_sym_parser();
   test_molp_method_parser();
   test_molp_wf_parser();
   test_molp_cas_reader();
   test_molp_basis_parser();
//   test_bessel_gaussian_poly_integral();
//   test_prim_radial_ovlp();
  // test_prim_ovlp();
   test_ao_ovlp();
   
   return 0;
}
