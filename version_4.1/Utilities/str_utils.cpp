#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include "utilities.h"

charTypeT charType(char c)
{
    if(isdigit(c))return digit;
    if(isalpha(c))return alpha;
    return other;
}

std::string separateThem(std::string inString)
{
   using namespace std;
  string oString = "";charTypeT st=other;
    for(unsigned int i=0;i!=inString.size();i++)
    {
       char c = inString[i];
        if( (st==alpha && charType(c)==digit) || (st==digit && charType(c)==alpha) )
          return oString;
         // oString.push_back(' ');
        oString.push_back(c);st=charType(c);
    }
    return 0;
}
int l_number(std::string bas_func_type)
{
   if(bas_func_type == "1s")
   {
      return 0;
   }
   else if(bas_func_type=="2px" || bas_func_type=="2py" || bas_func_type=="2pz" )
   {
      return 1;
   }
   else if(bas_func_type=="3d2-" || bas_func_type=="3d1-" || bas_func_type=="3d0" || bas_func_type=="3d1+" || bas_func_type=="3d2+" )
   {
      return 2;
   }
   else if( bas_func_type=="4f3-" || bas_func_type=="4f2-" || bas_func_type=="4f1-" || bas_func_type=="4f0" || bas_func_type=="4f1+" || bas_func_type=="4f2+" || bas_func_type=="4f3+")
   {
      return 3;
   }
   else if(bas_func_type=="5g4-" || bas_func_type=="5g3-" || bas_func_type=="5g2-" || bas_func_type=="5g1-" || bas_func_type=="5g0" || bas_func_type=="5g1+" || bas_func_type=="5g2+" || bas_func_type=="5g3+" || bas_func_type=="5g4+")
   {
      return 4;
   }
   else
   {
      std::cout<<"ERROR IN L VALUE DETERMINATION : "<<bas_func_type.c_str()<<std::endl;
      exit(EXIT_SUCCESS);
   }
}
int ml_number(std::string bas_func_type,int l)
{
   switch(l) {

   case 0:
      return 0;
   case 1:
      if(bas_func_type=="2px")
         return 1;
      else if(bas_func_type=="2pz")
         return 0;
      else if(bas_func_type=="2py")
         return -1;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   case 2:
      if(bas_func_type=="3d2-")
         return -2;
      else if(bas_func_type=="3d1-")
         return -1;
      else if(bas_func_type=="3d0")
         return 0;
      else if(bas_func_type=="3d1+")
         return 1;
      else if(bas_func_type=="3d2+")
         return 2;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   case 3:
      if(bas_func_type=="4f3-")
         return -3;
      else if(bas_func_type=="4f2-")
         return -2;
      else if(bas_func_type=="4f1-")
         return -1;
      else if(bas_func_type=="4f0")
         return 0;
      else if(bas_func_type=="4f1+")
         return 1;
      else if(bas_func_type=="4f2+")
         return 2;
      else if(bas_func_type=="4f3+")
         return 3;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   case 4:
      if(bas_func_type=="5g4-")
         return -4;
      if(bas_func_type=="5g3-")
         return -3;
      else if(bas_func_type=="5g2-")
         return -2;
      else if(bas_func_type=="5g1-")
         return -1;
      else if(bas_func_type=="5g0")
         return 0;
      else if(bas_func_type=="5g1+")
         return 1;
      else if(bas_func_type=="5g2+")
         return 2;
      else if(bas_func_type=="5g3+")
         return 3;
      else if(bas_func_type=="5g4+")
         return 4;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   default:
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);

   }
}

bool csf_string_parser(int n_sym_occ,int ci_size,std::vector<int> n_occ,std::vector<int> n_closed,std::vector<int> n_frozen, std::vector<std::string> csf_string,std::vector<int>* csf_mo,std::vector<int>* csf_spin)
{

   
   //This function parses the csf string that are gathered in molpro.
   //The strings come in a vector that has n_sym_occ*ci_size elements. In all the symetries, there is the same number of mo's and electrons


   // !!! Warning. We assume that the closed orbitals do not put a IR in the unoccupied category. here, n_sym_occ really means that the IR
   // orbitals are represented in the CSF in molpro

   std::string tmp_str;
   int mo_index;
   int elec_index;

   //We start by looping over the different CSF
   for (int i=0; i!=ci_size; i++)
   {
      //We represent the CSF by the state of its electrons. We have a spin part and a spatial part (mo)
      //The routine reads the occupation numbers and translates
      mo_index=0;
      elec_index=0;
               for(int l=0;l!=n_sym_occ;l++)
               {
                  //This selects one of the strings : csf_string[n_symocc*i+l]
                  tmp_str=csf_string[n_sym_occ*i+l];
                  if(n_closed[l]!=0)
                  {
                     for(int j=0;j!=n_closed[l]+n_frozen[l];j++)
                     {
                        csf_mo->push_back(mo_index);
                        csf_spin->push_back(0); // spin up
                        elec_index++;
                        csf_mo->push_back(mo_index);
                        csf_spin->push_back(1); // spin down
                        elec_index++;
                        mo_index++;
                     }
                 }
                 for(int j=0;j!=n_occ[l]-n_closed[l]-n_frozen[l];j++)
                 {
                    if(tmp_str.at(j)=='0')
                    {
                       mo_index++;
                       continue;
                    }
            
                    else if(tmp_str.at(j)=='2')
                    {
                        csf_mo->push_back(mo_index);
                        csf_spin->push_back(0); // spin up
                        elec_index++;
                        csf_mo->push_back(mo_index);
                        csf_spin->push_back(1); // spin down
                        elec_index++;
                        mo_index++;
                        continue;
                    }
            
                    else if(tmp_str.at(j)=='a' || tmp_str.at(j)=='/')
                    {
                        csf_mo->push_back(mo_index);
                        csf_spin->push_back(0); // spin up
                        elec_index++;
                        mo_index++;
                        continue;
                    }
             
                    else if( tmp_str.at(j) == 'b' || tmp_str.at(j) == '\\')
                    {
                       csf_mo->push_back(mo_index);
                       csf_spin->push_back(1); // spin down
                       elec_index++;
                       mo_index++;
                       continue;
                    }
                    else 
                       err_ci_not_rec_symbol(tmp_str,tmp_str.at(j));
                 }
               }

   }
      return 0;
}
