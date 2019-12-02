/*
std::complex<double> bessel_cont_angleint(int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,std::complex<double> **result)
{
   / *
    * For a given set of numbers l1,m1 and l2,m2, this routine will fill the *result array with the values
    * of the corresponding angular integrals. Namely, it will compute the angular part of the 
    *  - Bessel / Contraction overlap
    *  - Bessel / d/dx Contraction overlap
    *  - Bessel / d/dy Contraction overlap
    *  - Bessel / d/dz Contraction overlap
    * The cartesian components of this gradient are gven by angular integrals including the proper 
    * transformation of coordinates.
    * d/dx = sin(thet)*cos(phi)*d/dk + (1/k)*cos(thet)*cos(phi)*d/dthet - (1/(k*sin(thet)))*sin(phi)*d/dphi 
    * d/dy = sin(thet)*sin(phi)*d/dk + (1/k)*cos(thet)*sin(phi)*d/dthet + (1/(k*sin(thet)))*cos(phi)*d/dphi 
    * d/dz = cos(thet)*d/dk - (1/k)*sin(thet)*d/dthet 
    *
    * There are thus 9 integrals (including the overlap) to be computed, for every value of l3 [|l1-l2|,l1+l2] and m3 [-l3,l3]
    *
    * The structure of the result array is [n_int][n_l3][n_m3]
    *
    * [n_int] is such as follows : 
    * 0: overlap
    * 1: ddx - k component
    * 2 : ddx - thet component
    * 3: ddx - phi component
    * 4: ddy - k component
    * 5 : ddy - thet component
    * 6: ddy - phi component
    * 7: ddz - k component
    * 8 : ddz - thet component
    *
    * /
   int l1(jl);
   int l2(angular_mom_numbers[0]);
   int m1(jml);
   int m2(angular_mom_numbers[1]);

   
   for(int l3=int(fabs(l1-l2));l3<l1+l2+1;l3++)
   {
      for(int m3=-l3;m3<l3+1;m3++)
      {
         result[0][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*Dint(l1,l2,l3,m1,m2,m3);

         result[1][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_sph[0])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               * 0.5 * ( azim_integ(m1,m2,m3+1) + azim_integ(m1,m2,m3-1) )
               * (1./(2.*l3+1.)) * ( 
                     gaunt_formula(l1,l2,l3+1,fabs(m1),fabs(m2),fabs(m3)+1) 
                   - gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)+1) 
                 );
         result[4][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_sph[0]) * prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               * 0.5 * ( azim_integ(m1,m2,-m3-1) - azim_integ(m1,m2,-m3+1) )
               * (1./(2.*l3+1.)) * ( 
                     gaunt_formula(l1,l2,l3+1,fabs(m1),fabs(m2),fabs(m3)+1) 
                   - gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)+1) 
                 );
         result[7][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2]) * prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               * ( azim_integ(m1,m2,m3))
               * (-1./(2.*l3+1.)) * ( 
                     (l3-m3+1)*gaunt_formula(l1,l2,l3+1,fabs(m1),fabs(m2),fabs(m3)) 
                   + (l3+m3)*gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)) 
                 );
         if(m1 != 0 )
         {
            result[3][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *(m1/fabs(m1))
               *(0.5  * (-m2) * ( azim_integ(-m1+1,-m2,m3) - azim_integ(-m1-1,-m2,m3) ) 
                  +0.5 * (-m3) * (azim_integ(-m1+1,m2,-m3) - azim_integ(-m1-1,m2,-m3))
                )
               * (1./(2.*m1)) * ( 
                  gaunt_formula(l1-1,l2,l3,fabs(m1)+1,fabs(m2),fabs(m3)) 
                  + (l1+m1)*(l1+m1-1)*gaunt_formula(l1-1,l2,l3,fabs(m1)-1,fabs(m2),fabs(m3)) 
                                );    
            result[6][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *(0.5  * (-m2) * ( azim_integ(m1+1,-m2,m3) + azim_integ(m1-1,-m2,m3) ) 
                  +0.5 * (-m3) * (azim_integ(m1+1,m2,-m3) + azim_integ(m1-1,m2,-m3))
                )
               * (1./(2.*m1)) * ( 
                  gaunt_formula(l1-1,l2,l3,fabs(m1)+1,fabs(m2),fabs(m3)) 
                  + (l1+m1)*(l1+m1-1)*gaunt_formula(l1-1,l2,l3,fabs(m1)-1,fabs(m2),fabs(m3)) 
                                );    

         }
         else if ( m2 != 0 )
         {
            result[3][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *( (-m2) * ( azim_integ(-1,-m2,m3) ) + (-m3) * (azim_integ(-1,m2,-m3)) );
               * (1./(2.*m2)) * ( 
                  gaunt_formula(l1,l2-1,l3,fabs(m1),fabs(m2)+1,fabs(m3)) 
                  + (l2+m2)*(l2+m2-1)*gaunt_formula(l1,l2-1,l3,fabs(m1),fabs(m2)-1,fabs(m3)) 
                                );    
            result[6][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *( (-m2) * ( azim_integ(1,-m2,m3) ) + (-m3) * (azim_integ(-1,m2,-m3)) );
               * (1./(2.*m2)) * ( 
                  gaunt_formula(l1,l2-1,l3,fabs(m1),fabs(m2)+1,fabs(m3)) 
                  + (l2+m2)*(l2+m2-1)*gaunt_formula(l1,l2-1,l3,fabs(m1),fabs(m2)-1,fabs(m3)) 
                                );    
         }
         else if(m3 != 0)
         {
            result[3][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *( (-m2) * ( azim_integ(-1,-m2,m3) ) + (-m3) * (azim_integ(-1,m2,-m3)) );
               * (1./(2.*m3)) * ( 
                  gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)+1) 
                  + (l3+m3)*(l3+m3-1)*gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)-1) 
                                );    
            result[6][l3][m3]=rYlm(l3,m3,nucl_spher_pos[1],nucl_spher_pos[2])*prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3)
               *( (-m2) * ( azim_integ(1,-m2,m3) ) + (-m3) * (azim_integ(1,m2,-m3)) );
               * (1./(2.*m3)) * ( 
                  gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)+1) 
                  + (l3+m3)*(l3+m3-1)*gaunt_formula(l1,l2,l3-1,fabs(m1),fabs(m2),fabs(m3)-1) 
                                );    
         }
         else 
         {
            result[3][l3][m3]=0;
         }




   }
}
*/
double azim_integ(int m1,int m2,int m3)
{
   using namespace std;
   double sum(0);
   if(m1 ==0 && m2 ==0 && m3 ==0)
      return 2*acos(-1);
   else if(m1<0 && m2<0 && m3<0)
      return 0;
   else if(m1<0 && m2<0 && m3>=0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0)+bool(fabs(m1)-fabs(m2)+fabs(m3)==0)-bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else if(m1<0 && m2>=0 && m3<0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0)-bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else if(m1<0 && m2>=0 && m3>=0)
      return 0;
   else if(m1>=0 && m2<0 && m3<0)
      return 0.5*acos(-1)*(bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0)-bool(fabs(m1)-fabs(m2)-fabs(m3)==0)-bool(fabs(m1)+fabs(m2)+fabs(m3)==0));
   else if(m1>=0 && m2<0 && m3>=0)
      return 0;
   else if(m1>=0 && m2>=0 && m3<0)
      return 0;

   else if(m1 >=0 && m2>=0 && m3 >=0)
      return 0.5*acos(-1)*(bool(fabs(m1)+fabs(m2)+fabs(m3)==0)+bool(fabs(m1)-fabs(m2)-fabs(m3)==0)+bool(fabs(m1)-fabs(m2)+fabs(m3)==0)+bool(fabs(m1)+fabs(m2)-fabs(m3)==0));
   else 
   {
      std::cout<<"ERROR AZIM INTEGRAL TYPE CANNOT BE DEFINED "<<m1<<","<<m2<<","<<m3<<std::endl;
      exit(EXIT_SUCCESS);
   }
   /*
   if((m1+m2+m3) % 2 == 0)
   {

      if(m1 ==0 && m2 ==0 && m3 ==0)
         return 2*acos(-1);
      else if((m1<0 && m2<2 && m3<0 )|| (m1>0 && m2>0 && m3>0 ))
         return 0;
      else if(m1<0 && m2 < 0)
         sum=0.5*(bool(m3==(m1-m2))-bool(m3==(m1+m2)));
      else if(m1<0 && m2>=0)
         sum=0.5*(bool(m3==-fabs(m1-m2))+bool(m3==-fabs(m1+m2)));
      else if(m1>=0 && m2<0)
         sum=0.5*(bool(m3==-fabs(m1-m2))-bool(m3==-fabs(m1+m2)));
      else
         sum=0.5*(bool(m3==fabs(m1-m2))+bool(m3==fabs(m1+m2)));

      return 2*acos(-1)*sum;
   }
   else
      return 0;
      */
}
double rYlm (int l,int m,double thet,double phi)
{
//   std::cout<<"Entering rYlm with parameters "<<l<<","<<m<<","<<thet<<","<<phi<<std::endl;
      if(m<0)
      {
//      std::cout<<sqrt(2)*associated_legendre(l,-m,cos(thet))<<" * "<<sin(-m*phi)<<std::endl;
         return sqrt(2)*associated_legendre(l,-m,cos(thet))*sin(fabs(m)*phi);
      }
      else if(m>0)
      {
//      std::cout<<sqrt(2)*associated_legendre(l,m,cos(thet))<<" * "<<cos(-m*phi)<<std::endl;
         return sqrt(2)*associated_legendre(l,m,cos(thet))*cos(m*phi);
      }
      else
      {
//       std::cout<<associated_legendre(l,m,cos(thet)) << "*1"<<std::endl;
         return associated_legendre(l,m,cos(thet));
      }
}
double prefactor_rYlm(int l,int m)
{
//   double sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));
   double temp(1);
   double val(0);

   if(fabs(m) > fabs(l))
   {
      std::cout<<"FATAL ERROR IN SPHERICAL HARMONICS PREFACTOR. M>L:"<<m<<">"<<l<<std::endl;
      return 0;
//      exit(EXIT_SUCCESS);
   }

   if(m == 0)
   {
      return sqrt((2*l+1)/(4*acos(-1)));
   }
   else if(m > 0)
   {
      val= pow(-1,m) * sqrt(2.) * sqrt((2*l+1)/ (4*acos(-1)));
      temp=1.;
      for (int tt=l-m+1;tt!=l+m+1;tt++)
      {
         temp/=double(tt);
      }
      val*=sqrt(temp);
/*      if(val<0)
         std::cout<<"WOW ! THAT'S NOT NORMAL ! "<<std::endl;*/
      return val;
//      return sign * sqrt(2) * sqrt((2*l+1) * exp(0.5*(ln_factorial(l-m,lnfact_memo)-ln_factorial(l+m,lnfact_memo))) / (4*acos(-1)));
//      return sign * sqrt(2) * sqrt((2*l+1) * double(factorial(l-m,fact_memo)) / (4*acos(-1) * double(factorial(l+m,fact_memo))));
   }
   else 
   {
      /*temp=1.;
      for (int tt=l-fabs(m)+1;tt!=l+fabs(m)+1;tt++)
      {
         temp/=double(tt);
      }
      val*=temp;*/
      return prefactor_rYlm(l,fabs(m));
//      return sign * sqrt(2) * sqrt((2*l+1) * exp((ln_factorial(l+m,lnfact_memo)-ln_factorial(l-m,lnfact_memo)/2)) / (4*acos(-1)));
//      return sign * sqrt(2) * sqrt((2*l+1) * double(factorial(l+m,fact_memo)) / (4*acos(-1) * double( factorial(l-m,fact_memo))));
//      return sign * factorial(l+m)*prefactor_rYlm(l,-m)/factorial(l-m);
   }
}
double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3)
{
   double sum(0);

   if(l1==0 && l2 == 0 && l3 == 0 && m1 == 0 && m2 == 0 && m3 == 0 )
      return (acos(-1));

   else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;

   else if(l1>0)
   {
      if(l1 == m1)
         return -(2*l1-1)*gaunt_formula(l1-1,l2,l3,m1-1,m2,m3);

      else if(l1 == m1+1 && l1>=2)
         return -(2*m1+1)*gaunt_formula(m1,l2,l3,m1-1,m2,m3);

      else if(l1==1 && m1==0)
      {
         return (1./(2.*l2+1.))*((l2-m2+1)*J_int_m2(0,l2+1,l3,0,m2,m3)+(l2+m2)*J_int_m2(0,l2-1,l3,0,m2,m3));
      }
      if(l1>=m1+2)
         return (1./(double(l1-m1)*double(l1-m1-1)))*(double(2*l1-1)*gaunt_formula(l1-1,l2,l3,m1+1,m2,m3)
               +double(l1+m1)*(l1+m1-1)*J_int_m2(l1-2,l2,l3,m1,m2,m3));
/*
      if(l1-1>=m1+1)
         sum+=(double(2*l1-1)/double((l1-m1)*(l1-1-m1)))*gaunt_formula(l1-1,l2,l3,m1+1,m2,m3);
      if(l1-1>m1 && l1>=2)
         sum+=(double((2*l1-1)*(l1-1+m1))/double((l1-m1)*(l1-1-m1)))*J_int_m2(l1-2,l2,l3,m1,m2,m3);
      if(l1>=2)
         sum+=(double((l1-1+m1))/double((l1-m1)))*J_int_m2(l1-2,l2,l3,m1,m2,m3);
*/
   }
   else if(l2>0)
   {
      if(l2 == m2)
         return -(2*l2-1)*gaunt_formula(l1,l2-1,l3,m1,m2-1,m3);
      else if(l2 == m2+1 && l2>=2)
         return -(2*m2+1)*gaunt_formula(l1,m2,l3,m1,m2-1,m3);
      else if(l2==1 && m2==0)
      {
         return (1./(2.*l3+1.))*((l3-m3+1)*J_int_m2(0,0,l3+1,0,0,m3)+(l3+m3)*J_int_m2(0,0,l3-1,0,0,m3));
      }

      if(l2>=m2+2)
         return (1./(double(l2-m2)*double(l2-m2-1)))*(double(2*l2-1)*gaunt_formula(l1,l2-1,l3,m1,m2+1,m3)
               +double(l2+m2)*(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3));

      /*
      if(l2-1>=m2+1)
         sum+=(double(2*l2-1)/double((l2-m2)*(l2-1-m2)))*gaunt_formula(l1,l2-1,l3,m1,m2+1,m3);
      if(l2-1>m2 && l2>=2)
         sum+=(double((2*l2-1)*(l2-1+m2))/double((l2-m2)*(l2-1-m2)))*J_int_m2(l1,l2-2,l3,m1,m2,m3);
      if(l2>=2)
         sum+=(double((l2-1+m2))/double((l2-m2)))*J_int_m2(l1,l2-2,l3,m1,m2,m3);
*/
   }
   else
   {
      if(l3 == m3)
         return -(2*l3-1)*gaunt_formula(l1,l2,l3-1,m1,m2,m3-1);
      else if(l3 == m3+1 && l3>=2)
         return -(2*m3+1)*gaunt_formula(l1,l2,m3,m1,m2,m3-1);
      else if(l3==1 && m3==0)
         return 0;

      if(l3>=m3+2)
         return (1./(double(l3-m3)*double(l3-m3-1)))*(double(2*l3-1)*gaunt_formula(l1,l2,l3-1,m1,m2,m3+1)
               +double(l3+m3)*(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3));
/*
      if(l3-1>=m3+1)
         sum+=(double(2*l3-1)/double((l3-m3)*(l3-1-m3)))*gaunt_formula(l1,l2,l3-1,m1,m2,m3+1);
      if(l3-1>m3 && l3>=2)
         sum+=(double((2*l3-1)*(l3-1+m3))/double((l3-m3)*(l3-1-m3)))*J_int_m2(l1,l2,l3-2,m1,m2,m3);
      if(l3>=2)
         sum-=(double((l3-1+m3))/double((l3-m3)))*J_int_m2(l1,l2,l3-2,m1,m2,m3);
         */

   }
   std::cout<<"REACHED THE DEAD END OF J_INT_M2 FUNCTION. CASE IS"<<l1<<","<<m1<<";"<<l2<<","<<m2<<";"<<l3<<","<<m3<<std::endl;
   exit(EXIT_SUCCESS);
   return sum;
}
double J_int_m1(int l1,int l2,int l3,int m1,int m2,int m3)
{
   /*
   if((l1+m1+l2+m2+l3+m3)%2!=0) //If the integrand is odd
      return 0;
      */
//   else
/*   double norm1(1);
   double norm2(1);
   double temp(1);

   temp=1;
   for(int tt=fabs((l1+1)-(m1+1))+1;tt!=fabs((l1+1)+(m1+1)+1);tt++)
   {
      temp*=double(tt);
   }
   for(int tt=fabs(l1-m1)+1;tt!=fabs(l1+m1)+1;tt++)
   {
      temp/=double(tt);
   }
   norm2*=sqrt(temp);
   //std::cout<<"norm2 = "<<norm2<<std::endl;
*/
   if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;
   if(l1==0 || l1-1 < m1+1)
   {
      return -(1./(2.*l1+1))*gaunt_formula(l1+1,l2,l3,m1+1,m2,m3);
   }
   else
   {
      return (1./(2.*l1+1))*(gaunt_formula(l1-1,l2,l3,m1+1,m2,m3)-gaunt_formula(l1+1,l2,l3,m1+1,m2,m3));
   }
}
double J_int_p1(int l1,int l2,int l3,int m1,int m2,int m3)
{
   /*
   std::cout<<l1<<";"<<l2<<";"<<l3<<";"<<m1<<";"<<m2<<";"<<m3<<std::endl;
   std::cout<<gaunt_formula(l1+1,l2,l3,m1,m2,m3,fact_memo)<<std::endl;
   std::cout<<gaunt_formula(l1-1,l2,l3,m1,m2,m3,fact_memo)<<std::endl;
   std::cout<<"===="<<(l1-m1+1)*gaunt_formula(l1+1,l2,l3,m1,m2,m3,fact_memo)<<std::endl;
   std::cout<<"===="<<(l1+m1)*gaunt_formula(l1-1,l2,l3,m1,m2,m3,fact_memo)<<std::endl;
   std::cout<<"++++"<<(1./(2.*l1+1.))*((l1-m1+1)*gaunt_formula(l1+1,l2,l3,m1,m2,m3,fact_memo)-(l1+m1)*gaunt_formula(l1-1,l2,l3,m1,m2,m3,fact_memo))<<std::endl;
*/
//   if((l1+m1+l2+m2+l3+m3)%2==0) //If the integrand is odd
//      return 0;
//   else
   if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;
   if(l1==0 || l1-1 < m1)
   {
      return (1./(2.*l1+1.))*(double(l1-m1+1.)*gaunt_formula(l1+1,l2,l3,m1,m2,m3));
   }
   else
   {
      return (1./(2.*l1+1.))*(double(l1-m1+1.)*gaunt_formula(l1+1,l2,l3,m1,m2,m3)+double(l1+m1)*gaunt_formula(l1-1,l2,l3,m1,m2,m3));
   }
}
double J_int_D(int l1,int l2,int l3,int m1,int m2,int m3)
{
   if( l2 == 0 && l3 == 0  )
      return 0;

   else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;

   else if(l2 == 0)
      return (double(l3)/double(l3+m3+1))*(double(l3-m3+1)*J_int_m2(l1,l2,l3+1,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2,m3+1))-double(l3+m3)*J_int_m2(l1,l2,l3-1,m1,m2,m3);
   else if(l3 == 0)
      return (double(l2)/double(l2+m2+1))*(double(l2-m2+1)*J_int_m2(l1,l2+1,l3,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2+1,m3))-double(l2+m2)*J_int_m2(l1,l2-1,l3,m1,m2,m3);
   else
      return (double(l3)/double(l3+m3+1))*(double(l3-m3+1)*J_int_m2(l1,l2,l3+1,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2,m3+1))-double(l3+m3)*J_int_m2(l1,l2,l3-1,m1,m2,m3)
            +(double(l2)/double(l2+m2+1))*(double(l2-m2+1)*J_int_m2(l1,l2+1,l3,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2+1,m3))-double(l2+m2)*J_int_m2(l1,l2-1,l3,m1,m2,m3);
}
double J_int_m1_D(int l1,int l2,int l3,int m1,int m2,int m3)
{
//   if((l1+m1+l2+m2+l3+m3)%2==0) //If the integrand is odd
//      return 0;
//   else
   if( l2 == 0 && l3 == 0  )
      return 0;

   else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;

   else if(l2 == 0)
      return (1./(2.*l3+1.))*((l3+1)*(l3+m3)*gaunt_formula(l1,l2,l3-1,m1,m2,m3)-l3*(l3-m3+1)*gaunt_formula(l1,l2,l3+1,m1,m2,m3));
   else if(l3 == 0)
      return (1./(2.*l2+1.))*((l2+1)*(l2+m2)*gaunt_formula(l1,l2-1,l3,m1,m2,m3)-l2*(l2-m2+1)*gaunt_formula(l1,l2+1,l3,m1,m2,m3));
   else
      return (1./(2.*l2+1.))*((l2+1)*(l2+m2)*gaunt_formula(l1,l2-1,l3,m1,m2,m3)-l2*(l2-m2+1)*gaunt_formula(l1,l2+1,l3,m1,m2,m3))
            +(1./(2.*l3+1.))*((l3+1)*(l3+m3)*gaunt_formula(l1,l2,l3-1,m1,m2,m3)-l3*(l3-m3+1)*gaunt_formula(l1,l2,l3+1,m1,m2,m3));
/*
      return -(1./(4.*l1+2.))*(
         double(l3+m3)*double(l3-m3+1.)*gaunt_formula(l1-1,l2,l3,m1+1,m2,m3-1)
         -gaunt_formula(l1-1,l2,l3,m1+1,m2,m3+1)
         +double(l2+m2)*double(l2-m2+1)*gaunt_formula(l1-1,l2,l3,m1+1,m2-1,m3)
         -gaunt_formula(l1-1,l2,l3,m1+1,m2+1,m3)
         -double(l3+m3)*double(l3-m3+1)*gaunt_formula(l1+1,l2,l3,m1+1,m2,m3-1)
         +gaunt_formula(l1+1,l2,l3,m1+1,m2,m3+1)
         -double(l2+m2)*double(l2-m2+1)*gaunt_formula(l1+1,l2,l3,m1+1,m2-1,m3)
         +gaunt_formula(l1+1,l2,l3,m1+1,m2+1,m3)
         );
         */
}
double J_int_p1_D(int l1,int l2,int l3,int m1,int m2,int m3)
{
   if((l2 == 0 && l3 == 0))
      return 0;

   else if(l3 == 0)
   {
      return (double(l2+m2)/double(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3))
      -(double(l2)/double(2*l2+1))*(
            (double(l2-m2+1)/double(2*l2+3))*(double(l2-m2+2)*J_int_m2(l1,l2+2,l3,m1,m2,m3)+(l2+m2+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l2+m2)/(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3)));
      /*
         return (1/(2.*l1+1))*(
              (l1-m1+1)*(
                 (l2/(l2+m2+1))*( 
                    (l2-m2+1)*J_int_m2(l1+1,l2+1,l3,m1,m2,m3)
                    -gaunt_formula(l1+1,l2,l3,m1,m2+1,m3)) 
                 -(l2+m2)*J_int_m2(l1+1,l2-1,l3,m1,m2,m3))
             +(l1+m1)*(
                 (l2/(l2+m2+1))*( 
                    (l2-m2+1)*J_int_m2(l1-1,l2+1,l3,m1,m2,m3)
                    -gaunt_formula(l1-1,l2,l3,m1,m2+1,m3)) 
                 -(l2+m2)*J_int_m2(l1-1,l2-1,l3,m1,m2,m3)));
                 */
   }
   else if(l2 == 0)
   {
//      std::cout<<"l1+1 = "<<l1+1<<", l2 = "<<l2<<", l3+1 = "<<l3+1<<",m1 = "<<m1<<",m2 = "<<m2<<",m3 = "<<m3<<"====>"<<J_int_m2(l1+1,l2,l3+1,m1,m2,m3)<<std::endl;
      return (double(l3+m3)/double(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))
      -(double(l3)/double(2*l3+1))*(
            (double(l3-m3+1)/double(2*l3+3))*(double(l3-m3+2)*J_int_m2(l1,l2,l3+2,m1,m2,m3)+(l3+m3+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l3+m3)/(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3)));
/*         return (1./(2.*l1+1))*(
              (l1-m1+1)*(
                 (double(l3)/double(l3+m3+1))*double( 
                    (l3-m3+1)*J_int_m2(l1+1,l2,l3+1,m1,m2,m3)
                    -gaunt_formula(l1+1,l2,l3,m1,m2,m3+1)) 
                 -(l3+m3)*J_int_m2(l1+1,l2,l3-1,m1,m2,m3))
             +(l1+m1)*(
                 (double(l3)/double((l3+m3+1)))*( 
                    (l3-m3+1)*J_int_m2(l1-1,l2,l3+1,m1,m2,m3)
                    -gaunt_formula(l1-1,l2,l3,m1,m2,m3+1)) 
                 -(l3+m3)*J_int_m2(l1-1,l2,l3-1,m1,m2,m3)));
                 */
   }
   else
   {
      return (double(l2+m2)/double(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3))
      -(double(l2)/double(2*l2+1))*(
            (double(l2-m2+1)/double(2*l2+3))*(double(l2-m2+2)*J_int_m2(l1,l2+2,l3,m1,m2,m3)+(l2+m2+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l2+m2)/(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3)))
      + (double(l3+m3)/double(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))
      -(double(l3)/double(2*l3+1))*(
            (double(l3-m3+1)/double(2*l3+3))*(double(l3-m3+2)*J_int_m2(l1,l2,l3+2,m1,m2,m3)+(l3+m3+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l3+m3)/(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3)));
    /*     return (1./(2.*l1+1))*(
              (l1-m1+1)*(
                 (double(l3)/double(l3+m3+1))*( 
                    (l3-m3+1)*J_int_m2(l1+1,l2,l3+1,m1,m2,m3)
                    -gaunt_formula(l1+1,l2,l3,m1,m2,m3+1)) 
                 -(l3+m3)*J_int_m2(l1+1,l2,l3-1,m1,m2,m3))
             +(l1+m1)*(
                 (double(l3)/double(l3+m3+1))*( 
                    (l3-m3+1)*J_int_m2(l1-1,l2,l3+1,m1,m2,m3)
                    -gaunt_formula(l1-1,l2,l3,m1,m2,m3+1)) 
                 -(l3+m3)*J_int_m2(l1-1,l2,l3-1,m1,m2,m3))
              +(l1-m1+1)*(
                 (double(l2)/double(l2+m2+1))*( 
                    (l2-m2+1)*J_int_m2(l1+1,l2+1,l3,m1,m2,m3)
                    -gaunt_formula(l1+1,l2,l3,m1,m2+1,m3)) 
                 -(l2+m2)*J_int_m2(l1+1,l2-1,l3,m1,m2,m3))
             +(l1+m1)*(
                 (double(l2)/double(l2+m2+1))*( 
                    (l2-m2+1)*J_int_m2(l1-1,l2+1,l3,m1,m2,m3)
                    -gaunt_formula(l1-1,l2,l3,m1,m2+1,m3)) 
                 -(l2+m2)*J_int_m2(l1-1,l2-1,l3,m1,m2,m3)));
                 */
   }
}

double I_m1_integral(int m1,int m2,int m3)
{

   if(m1<0)
      return 0.5*(azim_integ(fabs(-m1-1),m2,m3)-azim_integ(fabs(-m1+1),m2,m3));
   else
   {
      if(m1 == 1)
         return 0.5*(azim_integ(-fabs(m1+1),m2,m3));
      else
         return 0.5*(azim_integ(-fabs(m1+1),m2,m3)+pow(-1,bool(m1-1<0))*azim_integ(-fabs(m1-1),m2,m3));
   }
}
double I_p1_integral(int m1,int m2,int m3)
{
   if(m1<0)
   {
      if(fabs(m1) == 1)
          return 0.5*(azim_integ(-2,m2,m3));
      else
          return 0.5*(azim_integ(-fabs(-m1+1),m2,m3)-pow(-1,bool(m1+1<0))*azim_integ(-fabs(m1+1),m2,m3));
   }
   else
      return 0.5*(azim_integ(fabs(m1+1),m2,m3)+azim_integ(fabs(-m1+1),m2,m3));

}
double I_m1_D_integral(int m1,int m2,int m3)
{
   if(m1<0)
         return -m2*0.5*(azim_integ(fabs(-m1-1),-m2,m3)-azim_integ(fabs(-m1+1),-m2,m3))
           -m3*0.5*(azim_integ(fabs(-m1-1),m2,-m3)-azim_integ(fabs(-m1+1),m2,-m3));
   else
   {
      if(m1 == 1)
         return -m2*0.5*(azim_integ(-fabs(m1+1),-m2,m3))-m3*0.5*((azim_integ(-fabs(m1+1),m2,-m3)));
      else
      {
         return -m2*0.5*(azim_integ(-fabs(m1+1),-m2,m3)-azim_integ(-fabs(m1-1),-m2,m3))
            -m3*0.5*((azim_integ(-fabs(m1+1),m2,-m3)-azim_integ(-fabs(m1-1),m2,-m3)));
      }
   }

}
double I_p1_D_integral(int m1,int m2,int m3)
{
   if(m1<0)
   {
      if(fabs(m1)==1)
         return -m2*0.5*(azim_integ(-fabs(-m1+1),-m2,m3))-m3*0.5*(azim_integ(-fabs(-m1+1),m2,-m3));
      else
      {
      return -m2*0.5*(-pow(-1,bool(m1+1<0))*azim_integ(-fabs(m1+1),-m2,m3)+azim_integ(-fabs(-m1+1),-m2,m3))
         -m3*0.5*(-pow(-1,bool(m1+1<0))*azim_integ(-fabs(m1+1),m2,-m3)+azim_integ(-fabs(-m1+1),m2,-m3));
      }
   }
   else
      return -m2*0.5*(azim_integ(fabs(m1+1),-m2,m3)+azim_integ(fabs(-m1+1),-m2,m3))
   -m3*0.5*(azim_integ(fabs(m1+1),m2,-m3)+azim_integ(fabs(-m1+1),m2,-m3));
}
double gaunt_formula(int l1,int l2,int l3,int m1,int m2,int m3)
{
//   if(l1 <= 12 && l2<=12 && l3<=12)
//      return GAUNT_INTEG_VAL[l1*l1+l1+m1][l2*l2+l2+m2][l3*l3+l3+m3];

//   if((l1+l2+l3)%2!=0)
//      return 0;
   double sign(1);
   double temp;
   double factor(0);
   double sum(0);
   double G12(0);
   double G123(0);
   double ALP_integ(0);

      if(m1<0)
      {
         /*
         temp=1;
         for (int tt=l1-fabs(m1)+1;tt!=l1+fabs(m1)+1;tt++)
         {
            temp/=double(tt);
         }
         sign=temp;//pow(-1.,m1)*temp;
         */
         m1=-m1;
      }
      if(m2<0)
      {
         /*
         temp=1;
         for (int tt=l2-fabs(m2)+1;tt!=l2+fabs(m2)+1;tt++)
         {
            temp/=double(tt);
         }
         sign*=temp;//pow(-1.,m2)*temp;
         */
         m2=-m2;
      }
      if(m3<0)
      {
         /*
         temp=1;
         for (int tt=l3-fabs(m3)+1;tt!=l3+fabs(m3)+1;tt++)
         {
            temp/=double(tt);
         }
         sign*=temp;//pow(-1.,m3)*temp;
         */
         m3=-m3;
      }

   int m12(m1+m2);
   int m123(m12+m3);

   int fac1[MAX_FACTORIAL_PRIME];
   double temp2(1);
   double val;

   if(l1 <0 || l2<0 || l3<0)
      return 0;

//      std::cout<<"Negative order in Gaunt integral"<<l1<<","<<l2<<","<<l3<<std::endl;
/*      if(m1 <0 || m2<0 || m3<0)
   {
      std::cout<<"Negative m value in Gaunt integral"<<m1<<","<<m2<<","<<m3<<std::endl;
      exit(EXIT_SUCCESS);
   }*/
   else if(m1>l1 || m2>l2 || m3>l3)
      return 0;
//   std::cout<<"GAUNT INTEGRAL  :  "<<l1<<","<<l2<<","<<l3<<" ; "<<m1<<","<<m2<<","<<m3<<std::endl;

      temp=1;
      for (int tt=l1-m1+1;tt!=l1+m1+1;tt++)
      {
         temp*=double(tt);
      }
      val=sqrt(temp);
//      std::cout<<"----"<<val;
      temp=1;
      for (int tt=l2-m2+1;tt!=l2+m2+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);
//      std::cout<<"*"<<sqrt(temp);
      temp=1;
      for (int tt=l3-m3+1;tt!=l3+m3+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

//      std::cout<<"*"<<sqrt(temp)<<" == ";
      double prefactor(val);

//      std::cout<<"prefactor : "<<prefactor<<std::endl;

//   double prefactor(exp(0.5*(ln_factorial(l1+m1,lnfact_memo)+ln_factorial(l2+m2,lnfact_memo)+ln_factorial(l3+m3,lnfact_memo)
//               -ln_factorial(l1-m1,lnfact_memo)-ln_factorial(l2-m2,lnfact_memo)-ln_factorial(l3-m3,lnfact_memo))));

   for(int l12=fabs(l1-l2);l12!=l1+l2+1;l12++)
   {
      if((l1+l2+l12)%2!=0 || m12 > l12)
         continue;
      else
      {

         G12=pow(-1.,m12)*(2.*l12+1.)*wigner3j(l1,l2,l12,0,0,0)*wigner3j(l1,l2,l12,m1,m2,-m12);
//         std::cout<<"l12 = "<<l12<<" ; 3J (l1 l2 l12 0 0 0) = "<<wigner3j(l1,l2,l12,0,0,0)<<std::endl;
//         std::cout<<"3J (l1 l2 l12 m1 m2 -m12 ) = "<<wigner3j(l1,l2,l12,m1,m2,-m12)<<" => G12 = "<<G12<<std::endl;
         for(int l123=fabs(l12-l3);l123!=l12+l3+1;l123++)
         {
            if((l12+l3+l123)%2!=0 || m123 > l123)
               continue;
            else
            {
               G123=pow(-1,m123)*(2.*l123+1.)*wigner3j(l12,l3,l123,0,0,0)*wigner3j(l12,l3,l123,m12,m3,-m123);
               //std::cout<<"l123 = "<<l123<<" -> G123 = "<<G123<<std::endl;
               temp=1;
               for (int tt=l123-m123+1;tt!=l123+m123+1;tt++)
               {
                  temp/=double(tt);
               }
               val=sqrt(temp);
               factor=val;
               //std::cout<<"temp 1 = "<<val<<std::endl;
//               factor=exp(0.5*(ln_factorial(l123-m123,lnfact_memo)-ln_factorial(l123+m123,lnfact_memo)));

               fact_prime_decomposer(int(double(l123-m123)/2),fac1);
               temp=1;
               for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
               {
                  temp*=pow(double(PRIME[i]),double(fac1[i]));
               }
//               std::cout<<int(double(l123-m123)/2)<<" => "<<temp<<std::endl;

               if(m123 == 0 && l123 == 0)
               {
                  ALP_integ=2.;
               }
               else
               {
                  ALP_integ=(pow(-1,m123)+pow(-1,l123))
                     *pow(2.,m123-2)*m123
                     *gamma_int_or_half(double(l123)/2.)
                     *gamma_int_or_half(double(l123+m123+1)/2.)
                     /(temp*gamma_int_or_half(double(l123+3)/2.));
//               std::cout<<l123<<","<<m123<<";"<<(pow(-1,m123)+pow(-1,l123))<<" * "<<pow(2.,m123-2)*m123<<" * "<<gamma_int_or_half(double(l123)/2.)<<"*"<<gamma_int_or_half(double(l123+m123+1)/2.)<<"/ ("<<gamma_int_or_half(double(l123+3)/2.)<<"*"<<temp<<")"<<" = "<<ALP_integ<<std::endl;
               }


//               std::cout<<G12<<"*"<<G123<<"*"<<factor<<"*"<<ALP_integ<<std::endl;
               sum+=G12*G123*factor*ALP_integ;
            }
         }
      }
   }
//   std::cout<<" = "<<prefactor*sum<<std::endl;
   return prefactor*sum;
}
/*double Dint(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   return sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*acos(-1)))*wigner3j(l1,l2,l3,0,0,0,lnfact_memo)*wigner3j(l1,l2,l3,m1,m2,m3,lnfact_memo);
}*/
