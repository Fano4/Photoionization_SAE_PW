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
   azim_integ_test();
   test_Jint_sort_indices();
   test_Jint_signflip_renormalize();
   test_Jint_normalize();
   test_Jint_special_cases();
   test_ALP_integral();
   return 0;
}
