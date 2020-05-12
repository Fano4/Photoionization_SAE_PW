#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

enum charTypeT{ other, alpha, digit};
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
