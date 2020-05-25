#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "filesutils.h"
#include <iomanip>
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
bool test_molp_basis_parser()
{
   bool test1(0);
   std::cout<<"Testing molp_basis_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";

   int tot_basis_size;
   std::vector<int> basis_size;
   std::vector<int> cont_num;
   std::vector<int> nuc_bas_func;
   std::vector<int> l_val;
   std::vector<int> m_val;
   std::vector<double> cont_zeta;
   std::vector<double> cont_coeff;

   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l_val,&m_val,&cont_zeta,&cont_coeff,input_file);

   std::cout<<" Basis size: "<<basis_size[0]<<","<<basis_size[1]<<","<<basis_size[2]<<","<<basis_size[3]<<std::endl;
   tot_basis_size=basis_size[0]+basis_size[1]+basis_size[2]+basis_size[3];
   int count(0);

   for(int i=0; i!=tot_basis_size;i++)
   {
      std::cout<<"Basis function "<<i<<" => "<<cont_num[i]<<" contractions. angular numbers ="<<l_val[i]<<","<<m_val[i]<<std::endl;
      for(int j=0;j!=cont_num[i];j++)
      {
         std::cout<<"coeff "<<cont_coeff[count]<<" ; zeta "<<cont_zeta[count]<<std::endl;
         count++;
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
}
bool test_molp_lcao_parser()
{
   using namespace std;
   bool test1(0);
   std::cout<<"Testing molp_lcao_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";
   int count(0);
   int n_sym;
   int method_index(0);
   std::vector<double> lcao_coeff;
   std::vector<int> n_occ;
   std::vector<int> n_closed;
   std::vector<int> n_frozen;
   std::vector<int> basis_size;
   std::vector<int> cont_num;
   std::vector<int> nuc_bas_func;
   std::vector<int> l_val;
   std::vector<int> m_val;
   std::vector<double> cont_zeta;
   std::vector<double> cont_coeff;


   molp_lcao_parser(method_index,&lcao_coeff,input_file);

   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file); 
   molp_basis_parser(&basis_size,&cont_num,&nuc_bas_func,&l_val,&m_val,&cont_zeta,&cont_coeff,input_file);
   n_sym=molp_sym_parser(input_file);

   std::cout<<"LCAO size : "<<lcao_coeff.size()<<std::endl;

   for(int s=0; s!=n_sym;s++)
   {
      std::cout<<"MO OF SYMMETRY "<<s+1<<std::endl<<std::endl;
      for(int mo=0;mo!=n_occ.at(s);mo++)
      {
         std::cout<<"MO "<<mo+1<<"."<<s+1<<" : "<<std::endl;
         for(int ao=0;ao!=basis_size[s];ao++)
         {
            if(ao%10 == 0 )
               cout<<endl;
            std::cout<<std::scientific<<std::setw(14)<<lcao_coeff.at(count);
            count++;
         }std::cout<<std::endl;
      }std::cout<<std::endl;
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
bool test_molp_ci_parser()
{
   using namespace std;
   bool test1(0);
   std::cout<<"Testing molp_ci_parser"<<std::endl;
   std::string input_file="/home/users/stephan/Photoionization_SAE_PW/version_4.1/Tests/LiH_6.325.out";
   int count(0);
   int countes(0);
   int n_sym;
   int method_index(0);
   std::vector<int> csf_mo;
   std::vector<int> csf_spin;
   std::vector<int> ci_num;
   std::vector<double> ci_coeff;
   std::vector<int> n_occ;
   std::vector<int> n_closed;
   std::vector<int> n_frozen;
   std::vector<int> n_elec;
   std::vector<int> sym;
   std::vector<int> spin;
   std::vector<int> charge;
   std::vector<int> n_states;


   molp_ci_parser(method_index,&csf_mo,&csf_spin,&ci_coeff,&ci_num,input_file);

   molp_cas_reader(method_index,&n_occ,&n_closed,&n_frozen,input_file); 
   n_sym=molp_sym_parser(input_file);
   molp_wf_parser(method_index,&n_elec,&sym,&spin,&charge,&n_states,input_file);

   std::cout<<"CI size : "<<ci_coeff.size()<<std::endl;

    std::cout<<"which means there are "<<std::endl;
   for(int s=0;s!=n_sym;s++)
   {
      std::cout<<ci_num.at(s)<<" CI vectors in symmetry "<<s+1<<std::endl;
   }

   for(int s=0; s!=n_sym;s++)
   {
      std::cout<<"CONFIGURATIONS OF SYMMETRY "<<s+1<<std::endl<<std::endl;
      //loop over the CIs in symmetry s
      for(int ci=0;ci<ci_num.at(s);ci++)
      {
         std::cout<<"CI "<<ci<<std::endl;
         //first display the configuration
         for(int m=0;m!=n_elec.at(s);m++)
            std::cout<<"( "<<csf_mo.at((ci+count)*n_elec.at(s)+m)<<" , "<<csf_spin.at((ci+count)*n_elec.at(s)+m)<<" ) ";
         for(int es=0;es!=n_states.at(s);es++)
         {
            if(es%8 == 0 && es != 0)
               cout<<endl;
            std::cout<<std::scientific<<std::setw(14)<<ci_coeff.at(countes+ci*n_states.at(s)+es);
         }std::cout<<std::endl;
      }std::cout<<std::endl;
      count+=ci_num.at(s);
      countes+=ci_num.at(s)*n_states.at(s);
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
