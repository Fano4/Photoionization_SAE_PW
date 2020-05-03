#include "photoion_header.h"

int main(int argc,char* argv[])
{
   if(argc<2)
   {
      // report version
      std::cout << argv[0] << " Version " << PhotoICE_VERSION_MAJOR << "."
              << PhotoICE_VERSION_MINOR << std::endl;
      std::cout << "Usage: " << argv[0] << " number" << std::endl;
      return 1;
   }
   //int photoion_comp(int argc, char* argv[]);
   //omp_set_num_threads(1); 
   int lmax=10;
   return 0;
}

