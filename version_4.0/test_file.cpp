#include "test_file.hpp"

void test_radial(int l1,int l2,int l3,int m1,int m2,int m3,double zeta,double kp,double* r0)
{ 
   //std::cout<<kp*kp*27.211/2<<","<<pow(-1,l1-l3)*4*acos(-1)*j_l(l3,kp*r0[0])*kp*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*rYlm(l3,m3,r0[1],r0[2])*Dint(l1,l2,l3,m1,m2,m3)<<std::endl;//<<","<<pow(0,-kp,l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*rYlm(l2,m2,r0[1],r0[2])<<std::endl;
//   std::cout<<kp*kp*27.211/2<<","<<std::abs(kp*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*j_l(l3,kp*r0))<<std::endl; 
//   std::cout<<kp*kp*27.211/2<<","<<std::abs((double(l2)/kp-kp/(2*zeta))*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2)*j_l(l3,kp*r0)+dj_ldz(l3,kp*r0)*pow(std::complex<double>(0,-kp),l2)*exp(-kp*kp/(4.*zeta))/pow(2.*zeta,1.5+l2))<<std::endl; 
   //std::cout<<j_l(l3,kp*r0)<<std::endl;
}
double test2_integral(int l1,int l2,int l3,int m1,int m2,int m3)
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
double analytic_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0)
{
   double Kval(pow(kp,l2)*exp(-kp*kp/(4*zeta))/pow(2*zeta,1.5+l2));
   double ddk_Kval((double(l2)/kp-kp/(2*zeta))*Kval);
   double bessel_val(0);
   double ddk_bessel_val(0);
   double ang_int7(0);
   double ang_int8(0);
   double ang_int9(0);
   double sum(0);

            for(int l3=0;l3!=l1+l2+2;l3++)//l3
            {
               for(int m3=-l3;m3!=l3+1;m3++)
               {
          //        std::cout<<J_int_p1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
          //        std::cout<<J_int_m1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
          //        std::cout<<azim_integ(m1,m2,m3)<<std::endl;
                  bessel_val=j_l(l3,kp*r0[0]);
   //               ddk_bessel_val=r0[0]*dj_ldz(l3,kp*r0[0]);
  /*                ang_int7=
                      pow(-1,((l2+l3-l1-1)/2))
                      *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2])
          //            *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                      *J_int_p1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*azim_integ(m1,m2,m3));

                  ang_int8=
                      pow(-1,((l2+l3-l1-1)/2))
                      *(-4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2])
        //              *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                      *J_int_m1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*azim_integ(m1,m2,m3));
                      */
                  ang_int9=
                         pow(-1,((l1-l2-l3)/2))
                          *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2]))
                          *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                          *gaunt_formula(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*azim_integ(m1,m2,m3);
//                  sum+=
//                     (kp*(ddk_Kval*bessel_val+Kval*ddk_bessel_val)*ang_int7
//                     +Kval*bessel_val*ang_int8)/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2*zeta,1.5+l2));

                  sum+=kp*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2*zeta,1.5+l2));
               }
            }
   std::cout<<"++"<<kp*kp*27.211/2<<"  "<<sum<<std::endl;
   return sum;
}
double numerical_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0)
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

   int nr(200);
   int nt(200);
   int nf(200);

   double R(10);
   double dr(R/nr);
   double dt(acos(-1)/nt);
   double df(2.*acos(-1)/nf);

   double bessel_norm=1;//kp*kp*pow(R,3)*(j_l(l1,kp*R)*j_l(l1,kp*R)-j_l(l1-1,kp*R)*j_l(l1+1,kp*R))/(acos(-1)));

   double x0(r0[0]*sin(r0[1])*cos(r0[2]));
   double y0(r0[0]*sin(r0[1])*sin(r0[2]));
   double z0(r0[0]*cos(r0[1]));

   double sum(0);

   for(int i=0;i!=nr;i++)
   {
      std::cout<<"r = "<<r<<std::endl;
      r=i*dr;
      for(int j=0;j!=nt;j++)
      {
         thet=j*dt;
         for(int l=0;l!=nf;l++)
         {
            phi=l*df;
            x=r*sin(thet)*cos(phi);
            y=r*sin(thet)*sin(phi);
            z=r*cos(thet);
            xp=x-x0;
            yp=y-y0;
            zp=z-z0;
            cart_to_spher(&xp,&yp,&zp,&rp,&thetp,&phip);

            // integral dxdydz * k * z * sqrt(2/pi) * j_l1(kr) * Y_l1m1(theta,phi) * Y_l2m2(thetap,phip) * |r-r0|^l2 * exp(-zeta * |r-r0|^2)
          //  std::cout<<z<<"  "<<dx*dy*dz*sqrt(2./acos(-1))*kp*z*rYlm(l1,m1,thet,phi)*j_l(l1,kp*r)*rYlm(l2,m2,thetp,phip)*pow(rp,l2)*exp(-zeta*pow(rp,2))<<std::endl;
//            sum+=dx*dy*dz*pow(rYlm(l2,m2,thetp,phip)*pow(rp,l2)*exp(-zeta*pow(rp,2))/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2*zeta,1.5+l2)),2);
//            sum+=r*r*dr*pow(sqrt(2./acos(-1))*kp*j_l(l1,kp*r)/sqrt(bessel_norm),2);
//            sum+=r*r*sin(thet)*dr*dt*df*pow(sqrt(2./acos(-1))*kp*rYlm(l1,m1,thet,phi)*j_l(l1,kp*r)/sqrt(bessel_norm),2);
            sum+=r*r*sin(thet)*dr*dt*df*sqrt(2./acos(-1))*kp*rYlm(l1,m1,thet,phi)*j_l(l1,kp*r)/sqrt(bessel_norm)*rYlm(l2,m2,thetp,phip)*pow(rp,l2)*exp(-zeta*pow(rp,2))/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2*zeta,1.5+l2));
            //sum+=dx*dy*dz*sqrt(2./acos(-1))*kp*z*rYlm(l1,m1,thet,phi)*j_l(l1,kp*r)*rYlm(l2,m2,thetp,phip)*pow(rp,l2)*exp(-zeta*pow(rp,2))/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2*zeta,1.5+l2));
         }
         //std::cout<<"===>"<<sum<<std::endl;
         //exit(EXIT_SUCCESS);
      }
   }
   std::cout<<kp*kp*27.211/2<<"  "<<sum<<std::endl;
}
void pw_bessel_comparison(double kp,double kthet,double kphi,double r,double thet,double phi)
{
   double kx(kp*sin(kthet)*cos(kphi));
   double ky(kp*sin(kthet)*sin(kphi));
   double kz(kp*cos(kthet));
   double x(r*sin(thet)*cos(phi));
   double y(r*sin(thet)*sin(phi));
   double z(r*cos(thet));

   double bessel_val(0);

    std::cout<<exp(std::complex<double>(0,kx*x+ky*y+kz*z))/pow(2*acos(-1),1.5)<<std::endl;

    std::complex<double> sum(0);
    std::complex<double> prefactor=0;
    double factor(0);

    for(int l1=0;l1!=15;l1++)
    {
       for(int m1=-l1;m1!=l1+1;m1++)
       {
          prefactor=pow(std::complex<double>(0,1),l1)*rYlm(l1,m1,kthet,kphi);
          factor=sqrt(2./acos(-1))*j_l(l1,kp*r)*rYlm(l1,m1,thet,phi);
          sum+=prefactor*factor;
       }
       std::cout<<"+++"<<sum<<std::endl;
    }
}
void pw_bessel_overlap_comparison(int l2,int m2,double zeta,double kp,double thet,double phi,double* r0)
{
   double kx(kp*sin(thet)*cos(phi));
   double ky(kp*sin(thet)*sin(phi));
   double kz(kp*cos(thet));

   double xp(r0[0]*sin(r0[1])*cos(r0[2]));
   double yp(r0[0]*sin(r0[1])*sin(r0[2]));
   double zp(r0[0]*cos(r0[1]));

   double Kval(pow(kp,l2)*exp(-kp*kp/(4*zeta))/(pow(2*zeta,1.5+l2)));
   double bessel_val(0);
   std::complex<double> ang_int9(0);

    std::cout<<pow(std::complex<double>(0,-1),l2)*Kval*rYlm(l2,m2,thet,phi)*exp(std::complex<double>(0,-(kx*xp+ky*yp+kz*zp)))/sqrt(0.5*intplushalf_gamma(l2+1)/(pow(2.*zeta,1.5+l2)))<<std::endl;

    std::complex<double> sum(0);
    std::complex<double> prefactor=0;
    double factor(0);

    for(int l1=0;l1!=10;l1++)
    {
       for(int m1=-l1;m1!=l1+1;m1++)
       {
//          std::cout<<l1<<" , "<<m1<<" - "<<thet<<" , "<<phi<<" rYlm = "<<rYlm(l1,m1,thet,phi)<<std::endl;
          prefactor=pow(std::complex<double>(0,-1.),l1)*rYlm(l1,m1,thet,phi);
          for(int l3=fabs(l1-l2);l3!=l1+l2+1;l3++)
          {
             bessel_val=j_l(l3,kp*r0[0]);
//             std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" : j_"<<l3<<"("<<kp*r0[0]<<") = "<<bessel_val<<std::endl;
             for(int m3=-l3;m3!=l3+1;m3++)
             {
                  ang_int9=
                         pow(std::complex<double>(0,-1),l1-l2-l3)
                          *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2]))
                          *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                          *gaunt_formula(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*azim_integ(m1,m2,m3);

               /*   if(prefactor*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2)) != std::complex<double>(0,0))
                  {
                     std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" Sign = "<<pow(-1,(int(l1-l2+l3)/2))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Gaunt = "<<gaunt_formula(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Prefactor = "<<prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : rYl3m3 = "<<rYlm(l3,m3,r0[1],r0[2])<<std::endl;
                      std::cout<<"azim = "<<azim_integ(m1,m2,m3)<< " ; Bessel = "<< bessel_val<<std::endl;
                     std::cout<<" ===== "<<prefactor*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2))<<std::endl;
                  }*/
                  sum+=prefactor*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2));
             }
          }
          
       }
       std::cout<<"+++"<<sum<<std::endl;
    }
}
void pw_bessel_gradient_y_comparison(int l2,int m2,double zeta,double kp,double thet,double phi,double* r0)
{
   double kx(kp*sin(thet)*cos(phi));
   double ky(kp*sin(thet)*sin(phi));
   double kz(kp*cos(thet));

   double xp(r0[0]*sin(r0[1])*cos(r0[2]));
   double yp(r0[0]*sin(r0[1])*sin(r0[2]));
   double zp(r0[0]*cos(r0[1]));

   double Kval(pow(kp,l2)*exp(-kp*kp/(4.*zeta))/(pow(2*zeta,1.5+l2)));
   double ddk_Kval((double(l2)/kp-kp/(2.*zeta))*Kval);
   double bessel_val(0);
   double ddk_bessel_val(0);
   std::complex<double> ang_int4(0);
   std::complex<double> ang_int5(0);
   std::complex<double> ang_int6(0);
   std::complex<double> pw_integf(0);
   std::complex<double> f(sin(thet)*sin(r0[1])*cos(phi-r0[2])+cos(thet)*cos(r0[1]));
   std::complex<double> dfdt(cos(thet)*sin(r0[1])*cos(phi-r0[2])-sin(thet)*cos(r0[1])); 
   std::complex<double> dfdf(-sin(thet)*sin(r0[1])*sin(phi-r0[2]));
   std::complex<double> dfdf2(-sin(r0[1])*sin(phi-r0[2]));
   std::complex<double> phase_factor(exp(-std::complex<double>(0,1)*kp*r0[0]*f));
   std::complex<double> ddk_exp(-std::complex<double>(0,1)*r0[0]*f*phase_factor );
   std::complex<double> ddt_exp(-std::complex<double>(0,1)*r0[0]*dfdt*phase_factor);

    //std::cout<<pow(std::complex<double>(0,-1),l2-1)*Kval*rYlm(l2,m2,thet,phi)*exp(std::complex<double>(0,-(kx*xp+ky*yp+kz*zp)))/sqrt(0.5*intplushalf_gamma(l2+1)/(pow(2.*zeta,1.5+l2)))<<std::endl;

   std::complex<double> pw_integk(pow(std::complex<double>(0,-1),l2-1)*rYlm(l2,m2,thet,phi)*(ddk_Kval*phase_factor+Kval*ddk_exp)); 
   std::complex<double> pw_integt(pow(std::complex<double>(0,-1),l2-1) 
               *(ddt_exp*rYlm(l2,m2,thet,phi)+phase_factor*pow(2.,0.5*bool(m2!=0))*associated_legendre_der(l2,fabs(m2),cos(thet))*(bool(m2>=0) * cos(m2 * phi) + bool(m2<0) * sin(fabs(m2) * phi)))
               *Kval);

//   std::cout<<"***"<<sin(thet)*sin(phi)*pw_integk/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2))<<std::endl;
   std::cout<<"***"<<cos(thet)*sin(phi)*pw_integt/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2))<<std::endl;

   if(m2!=0)
   {
      pw_integf=(
         pow(2.,0.5*bool(m2!=0))*
         ((1./(2.*fabs(m2)))*sqrt(((2.*l2+1.)*factorial(l2-fabs(m2)))/((4.*acos(-1))*(factorial(l2+fabs(m2)))))
         *(associated_legendre_nonorm(l2+1,fabs(m2)+1,cos(thet)) + (l2-fabs(m2)+1) * (l2-fabs(m2)+2) * associated_legendre_nonorm(l2+1,fabs(m2)-1,cos(thet)))
         *(fabs(m2) * bool(m2>=0) * sin(m2 * phi) - fabs(m2) * bool(m2<0) * cos(m2 * phi)
         +std::complex<double>(0,1) * kp * r0[0] * dfdf // Real spherical harmonics 
         *(bool(m2>=0) * cos(m2* phi) + bool(m2<0) * sin(fabs(m2) * phi)))
          ) // ml sign dependent phi part
         *(pow(std::complex<double>(0,-1),l2-1) * Kval ) // Radial part
         *phase_factor); // phase factor
   }
   else
   {
      pw_integf=(
            std::complex<double>(0,1) * kp * r0[0] * associated_legendre(l2,m2,cos(thet))
            *(pow(std::complex<double>(0,-1),l2) * Kval)
            *phase_factor);
   }

//    std::cout<<pw_integf<<std::endl;
//    std::cout<<"=="<<(sin(thet)*sin(phi)*pw_integk+cos(thet)*sin(phi)*pw_integt+cos(phi)*pw_integf)/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2))<<std::endl;
          
    std::complex<double> sum(0);
    std::complex<double> prefactor=0;
    double factor(0);

    for(int l1=0;l1!=10;l1++)
    {
       for(int m1=-l1;m1!=l1+1;m1++)
       {
//          std::cout<<l1<<" , "<<m1<<" - "<<thet<<" , "<<phi<<" rYlm = "<<rYlm(l1,m1,thet,phi)<<std::endl;
          prefactor=pow(std::complex<double>(0,1.),1)*rYlm(l1,m1,thet,phi);
          for(int l3=fabs(l1-l2)-1;l3!=l1+l2+2;l3++)
          {
             if(l3<0)
                continue;
             bessel_val=j_l(l3,kp*r0[0]);
             ddk_bessel_val=r0[0]*dj_ldz(l3,kp*r0[0]);
//             std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" : j_"<<l3<<"("<<kp*r0[0]<<") = "<<bessel_val<<std::endl;
             for(int m3=-l3;m3!=l3+1;m3++)
             {
//                ang_int4=gaunt_formula(l1,l2,l3,m1,m2,m3);

//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<"prefactors : "<<prefactor_rYlm(l1,fabs(m1))<<" , "<<prefactor_rYlm(l2,fabs(m2))<<" , "<<prefactor_rYlm(l3,fabs(m3))<<std::endl;
//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Gaunt = "<<ang_int4<<std::endl;
//                  std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Jm1 = "<<J_int_m1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Jm2 = "<<J_int_m2(l1,l2,l3,m1,m2,m3)<<std::endl;
//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Jp1 = "<<J_int_p1(l1,l2,l3,m1,m2,m3)<<std::endl;
//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Jm1D = "<<J_int_m1_D(l1,l2,l3,m1,m2,m3)<<std::endl;
                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Jp1D = "<<J_int_p1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
                  std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Im1 = "<<I_m1_integral(m1,m2,m3)<<std::endl;
//                  std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Ipm1 = "<<I_p1_integral(m1,m2,m3)<<std::endl;
/*                ang_int4=
                         pow(std::complex<double>(0,-1),(l2+l3))
                          *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2]))
                          *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                          *J_int_m1(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*I_m1_integral(m1,m2,m3);
*/
//                std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : ang int 4 = "<<ang_int4<<std::endl;

                ang_int5=
                         pow(std::complex<double>(0,-1),(l2+l3))
                          *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2]))
                          *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                          *J_int_p1_D(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*I_m1_integral(m1,m2,m3);
 /*               ang_int6=
                          pow(std::complex<double>(0,-1),(l2+l3))
                          *(4.*acos(-1)*rYlm(l3,m3,r0[1],r0[2]))
                          *prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))
                          *J_int_m2(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))*I_p1_D_integral(m1,m2,m3);
*/
               /*   if(prefactor*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2)) != std::complex<double>(0,0))
                  {
                     std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" Sign = "<<pow(-1,(int(l1-l2+l3)/2))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Gaunt = "<<gaunt_formula(l1,l2,l3,fabs(m1),fabs(m2),fabs(m3))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : Prefactor = "<<prefactor_rYlm(l1,fabs(m1))*prefactor_rYlm(l2,fabs(m2))*prefactor_rYlm(l3,fabs(m3))<<std::endl;
                      std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : rYl3m3 = "<<rYlm(l3,m3,r0[1],r0[2])<<std::endl;
                      std::cout<<"azim = "<<azim_integ(m1,m2,m3)<< " ; Bessel = "<< bessel_val<<std::endl;
                     std::cout<<" ===== "<<prefactor*Kval*bessel_val*ang_int9/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2))<<std::endl;
                  }*/
                  sum+=prefactor*((1./kp)*Kval*bessel_val*ang_int5+(1./kp)*Kval*bessel_val*ang_int6)/sqrt(0.5*intplushalf_gamma(l2+1)/pow(2.*zeta,1.5+l2));
//                  sum+=prefactor*((ddk_Kval*bessel_val+Kval*ddk_bessel_val)*ang_int4);
             }
          }
          
       }
       std::cout<<"+++"<<sum<<std::endl;
    }
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
