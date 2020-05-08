#include "photoion_header.h"

#include "mathfunctions.h"

#include<cmath>

int main(int argc,char* argv[])
{
   if(argc<2)
   {
      // report version
      std::cout << argv[0] << " Version " << PhotoICE_VERSION_MAJOR << "."
              << PhotoICE_VERSION_MINOR << std::endl;

      return 1;
   }
   //int photoion_comp(int argc, char* argv[]);
   //omp_set_num_threads(1); 
   return 0;
}

