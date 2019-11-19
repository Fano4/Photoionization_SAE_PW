#include "test_file.hpp"

void test_radial(int l1,int l2,int l3,int m1,int m2,int m3,double zeta,double kp,double* r0,double* lnfact_memo)
{ 
   std::cout<<kp*kp*27.211/2<<","<<pow(-1,l1-l3)*4*acos(-1)*j_l(l3,kp*r0[0],lnfact_memo)*kp*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*rYlm(l3,m3,r0[1],r0[2],lnfact_memo)*Dint(l1,l2,l3,m1,m2,m3,lnfact_memo)<<std::endl;//<<","<<pow(0,-kp,l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*rYlm(l2,m2,r0[1],r0[2],lnfact_memo)<<std::endl;
//   std::cout<<kp*kp*27.211/2<<","<<std::abs(kp*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*j_l(l3,kp*r0,lnfact_memo))<<std::endl; 
//   std::cout<<kp*kp*27.211/2<<","<<std::abs((double(l2)/kp-kp/(2*zeta))*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*j_l(l3,kp*r0,lnfact_memo)+dj_ldz(l3,kp*r0,lnfact_memo)*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2))<<std::endl; 
   //std::cout<<j_l(l3,kp*r0,lnfact_memo)<<std::endl;
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
double analytic_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0,double* lnfact_memo)
{
   double Kval(pow(-kp,l2)*exp(-kp*kp/(4*zeta))/pow(2*zeta,1.5+l2));
   double ddk_Kval((double(l2)/kp-kp/(2*zeta))*Kval);
   double bessel_val(0);
   double ddk_bessel_val(0);
   double ang_int7(0);
   double ang_int8(0);
   double sum(0);

            for(int l3=0;l3!=l1+l2+2;l3++)//l3
            {
               for(int m3=-l3;m3!=l3+1;m3++)
               {
          //        std::cout<<J_int_p1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3),lnfact_memo)<<std::endl;
          //        std::cout<<J_int_m1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3),lnfact_memo)<<std::endl;
          //        std::cout<<azim_integ(m1,m2,m3)<<std::endl;
                  bessel_val=j_l(l3,kp*r0[0],lnfact_memo);
                  ddk_bessel_val=r0[0]*dj_ldz(l3,kp*r0[0],lnfact_memo);
                  ang_int7=
                      4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2],lnfact_memo)
                      *prefactor_rYlm(l1,fabs(m1),lnfact_memo)*prefactor_rYlm(l2,fabs(m2),lnfact_memo)*prefactor_rYlm(l3,fabs(m3),lnfact_memo)
                      *J_int_p1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3),lnfact_memo)*azim_integ(m1,m2,m3);

                  ang_int8=
                      -4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2],lnfact_memo)
                      *prefactor_rYlm(l1,fabs(m1),lnfact_memo)*prefactor_rYlm(l2,fabs(m2),lnfact_memo)*prefactor_rYlm(l3,fabs(m3),lnfact_memo)
                      *J_int_m1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3),lnfact_memo)*azim_integ(m1,m2,m3);
                  sum+=
                     (kp*(ddk_Kval*bessel_val+Kval*ddk_bessel_val)*ang_int7
                     +Kval*bessel_val*ang_int8)/sqrt(0.5*intplushalf_gamma(l2+1,lnfact_memo)/pow(2*zeta,1.5+l2));
               }
            }
   std::cout<<"++"<<kp*kp*27.211/2<<"  "<<sum<<std::endl;
   return sum;
}
double numerical_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0,double* lnfact_memo)
{
   //Computes the numerical value of the transition dipole between an atomic orbital centered at r0 with exponential factor zeta and a basis function sqrt(2*k^2/pi)*j_l(kr)*Y_lm(theta,phi) IN CARTESIAN COORDINATES
   // integral dxdydz * k * z * sqrt(2/pi) * j_l1(kr) * Y_l1m1(theta,phi) * Y_l2m2(thetap,phip) * |r-r0|^l2 * exp(-zeta * |r-r0|^2)

   double x(0);
   double y(0);
   double z(0);
   double xp(0);
   double yp(0);
   double zp(0);

   double r(0);
   double rp(0);
   double thet(0);
   double thetp(0);
   double phi(0);
   double phip(0);

   int nx(200);
   int ny(200);
   int nz(400);

   double xmin(-10);
   double ymin(-10);
   double zmin(-10);
   double xmax(10);
   double ymax(10);
   double zmax(10);

   double dx((xmax-xmin)/nx);
   double dy((ymax-ymin)/ny);
   double dz((zmax-zmin)/nz);

   double x0(r0[0]*sin(r0[1])*cos(r0[2]));
   double y0(r0[0]*sin(r0[1])*sin(r0[2]));
   double z0(r0[0]*cos(r0[1]));

   double sum(0);

   for(int i=0;i!=nx+1;i++)
   {
      x=xmin+i*(xmax-xmin)/nx;
      xp=x-x0;
      for(int j=0;j!=ny+1;j++)
      {
         y=ymin+j*(ymax-ymin)/ny;
         yp=y-y0;
         for(int l=0;l!=nz+1;l++)
         {
            z=zmin+l*(zmax-zmin)/nz;
            zp=z-z0;
            cart_to_spher(&x,&y,&z,&r,&thet,&phi);
            cart_to_spher(&xp,&yp,&zp,&rp,&thetp,&phip);

            // integral dxdydz * k * z * sqrt(2/pi) * j_l1(kr) * Y_l1m1(theta,phi) * Y_l2m2(thetap,phip) * |r-r0|^l2 * exp(-zeta * |r-r0|^2)
          //  std::cout<<z<<"  "<<dx*dy*dz*sqrt(2./acos(-1))*kp*z*rYlm(l1,m1,thet,phi,lnfact_memo)*j_l(l1,kp*r,lnfact_memo)*rYlm(l2,m2,thetp,phip,lnfact_memo)*pow(rp,l2)*exp(-zeta*pow(rp,2))<<std::endl;
            sum+=dx*dy*dz*sqrt(2./acos(-1))*kp*z*rYlm(l1,m1,thet,phi,lnfact_memo)*j_l(l1,kp*r,lnfact_memo)*rYlm(l2,m2,thetp,phip,lnfact_memo)*pow(rp,l2)*exp(-zeta*pow(rp,2))/sqrt(0.5*intplushalf_gamma(l2+1,lnfact_memo)/pow(2*zeta,1.5+l2));
         }
         //std::cout<<"===>"<<sum<<std::endl;
         //exit(EXIT_SUCCESS);
      }
   }
   std::cout<<kp*kp*27.211/2<<"  "<<sum<<std::endl;
}
void cart_to_spher(double* x,double* y,double* z,double * r,double* t,double *f)
{  
      *r=sqrt(*x * *x + *y * *y + *z * *z);
      
      if(*r==0)
      {  
         *t=0;
         *f=0;
      }
      else
      {  
         *t=acos(*z / *r);
         if(*x == 0 && *y > 0)
         {  
            *f=acos(-1)/2.;
         }
         else if (*x == 0 && *y < 0 )
         {  
            *f =3.*acos(-1)/2.;
         }
         else
         {  
            *f = atan2(*y,*x);
         }
      }
      if(*f < 0)
         *f+=2*acos(-1);
}
