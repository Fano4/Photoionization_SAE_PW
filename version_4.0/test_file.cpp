#include "test_file.hpp"

void test_radial(int l2,int l3,double zeta,double kp,double r0,double* lnfact_memo)
{
   std::cout<<pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)<<"   "; 
   std::cout<<j_l(l3,kp*r0,lnfact_memo)<<std::endl;
}
double test2_integral(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   double xmin(-1);
   double xmax(1);
   int nx(10000000);
   double dx((xmax-xmin)/double(nx));

   double x(0);
   double sum(0);

   for(int i=0;i!=nx+1;i++)
   {
      x=xmin+i*dx;
      sum+=pow(-1,m1+m2+m3)*associated_legendre_nonorm(l1,m1,x)*associated_legendre_nonorm(l2,m2,x)*associated_legendre_nonorm(l3,m3,x)*dx;
//      std::cout<<x<<" ; "<<associated_legendre_nonorm(l1,m1,x)*associated_legendre_nonorm(l2,m2,x)*associated_legendre_nonorm(l3,m3,x)<<std::endl;
   }

//   std::cout<<"===="<<sum<<std::endl;
   return sum;
}

