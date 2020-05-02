#include "test_headers.hpp"

//TEST FUNCTION PATTERN
/*
bool test_function()
{
   bool test1(0);

   std::cout<<" Testing function...";
   //Test_case_1
   {
     //Environment
     //test_case
     std::cout<<"1";
     //test
     if(test)
         test1=1;
      //Finalize
   }
   //Print results 
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }

}

*/

int main(int argc,char* argv[])
{
   test_determinant();
   test_legendre();
   test_prime_decomposer();
   test_fact_prime_decomposer();
   //test_factorized_sum();
   //test_wigner3j();
   test_associated_legendre_nonorm();
   gaunt_formula_test();
//   J_int_m2_test();
//   azim_integ_test();
   return 0;
}