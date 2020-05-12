#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "filesutils.h"
bool test_search()
{
   std::vector<int> match_pos;
   std::vector<int> num_of_match;
   bool found(0);

   found=search(&match_pos,&num_of_match,"wf","/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/search_test_files/test1.txt");
   if(found)
   {
      std::cout<<"Pattern found, "<<num_of_match[0]<<" occurences"<<std::endl;

      for(int i=0;i!=num_of_match[0];i++)
         std::cout<<"Position "<<match_pos[i]<<std::endl;
   }
   return 1;
}
bool test_molp_sym_parser()
{
   bool test1(0);
   std::cout<<"Testing molp_sym_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   int test_sym;
   int sym(4);

   test_sym=molp_sym_parser(input_file);

   if(sym==test_sym)
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
         std::cout<<"Error 1..."<<test_sym;
      }
      std::cout<<std::endl;
      return 0;
   }
}
bool test_molp_method_parser()
{
   bool test1(0);
   std::cout<<"Testing molp_method_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   int method(2);
   std::vector<int> test_method;
   molp_method_parser(&test_method,input_file);

   std::cout<<"Got method "<<test_method[0]<<" at position "<<test_method[1]<<std::endl;
   if(method==test_method[0])
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
         std::cout<<"Error 1..."<<test_method[0];
      }
      std::cout<<std::endl;
      return 0;
   }
}
bool test_molp_wf_parser()
{
   bool test1(0);
   std::cout<<"Testing molp_wf_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   int method_index(0);
   std::vector<int> n_elec;
   std::vector<int> sym;
   std::vector<int> spin;
   std::vector<int> charge;
   std::vector<int> n_states;

   molp_wf_parser(method_index,&n_elec,&sym,&spin,&charge,&n_states,input_file);

   for(unsigned int i=0;i!=n_elec.size();i++)
   {
      std::cout<<"Got wf "<<n_elec[i]<<","<<sym[i]<<","<<spin[i]<<";state,"<<n_states[i]<<std::endl;
      test1=1;
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
bool test_molp_cas_reader()
{
   bool test1(0);
   std::cout<<"Testing molp_cas_reader"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   int method_index(0);
   std::vector<int> n_occ;
   std::vector<int> n_closed;
   std::vector<int> n_frozen;

   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file);

   for(int i=0;i!=4;i++)
   {
      std::cout<<"Got cas "<<n_occ[i]<<","<<n_closed[i]<<std::endl;
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
}
