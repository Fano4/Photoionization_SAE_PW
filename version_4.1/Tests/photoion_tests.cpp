#include "photoion_tests.h"
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
   test_J_int_m1_D();
   test_J_int_p1_D();
   return 0;
}
