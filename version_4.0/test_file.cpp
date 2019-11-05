#include "test_file.hpp"

double test2_integral(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   double xmin(-1);
   double xmax(1);
   int nx(100000000);
   double dx((xmax-xmin)/nx);

   double x(0);
   double sum(0);

   for(int i=0;i!=nx+1;i++)
   {
      x=xmin+i*dx;
      sum+=associated_legendre_nonorm(l1,m1,x)*associated_legendre_nonorm(l2,m2,x)*associated_legendre_nonorm(l3,m3,x)*dx;
//      std::cout<<sum<<std::endl;
   }

//   std::cout<<"===="<<sum<<std::endl;
   return sum;
}

