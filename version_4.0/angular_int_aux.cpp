double azim_integ(int m1,int m2,int m3)
{
   //CHECKED ON FEB 28 2020
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
double polar_integ(int l1,int l2,int l3,int m1,int m2,int m3)
{
   // This function uses the formula for an integral over three spherical harmonics ant the azim_integ function to compute the integral over three legendre polynomials
   // S = A.J

   double prefactor(sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*acos(-1))));
   //double I(azim_integ(m1,m2,m3));
   double cg1(WignerSymbols::wigner3j(l1,l2,l3,0,0,0));
   double cg2(WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3));
   double A(prefactor_rYlm(l1,m1)*prefactor_rYlm(l2,m2)*prefactor_rYlm(l3,m3));

   return prefactor*cg1*cg2/(A);
}
double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3)
{
   // CHECKED ON MARCH 2 2020
   double sum(0);
   double sign(1);

   if(m1<0)
   {
      sign*=pow(-1,m1);
      m1=-m1;
   }
   if(m2<0)
   {
      sign*=pow(-1,m2);
      m2=-m2;
   }
   if(m3<0)
   {
      sign*=pow(-1,m3);
      m3=-m3;
   }

   if(l1==0 && l2 == 0 && l3 == 0 && m1 == 0 && m2 == 0 && m3 == 0 )
      return (acos(-1));

   else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;

   else if(l1>0)
   {
      if(l1 == m1)
         return -sign*((2*l1-1)*gaunt_formula(l1-1,l2,l3,m1-1,m2,m3));

      else if(l1 == m1+1 && l1>=2)
         return -sign*((2*m1+1)*gaunt_formula(m1,l2,l3,m1-1,m2,m3));

      else if(l1==1 && m1==0)
      {
         return sign*((1./(2.*l2+1.))*((l2-m2+1)*J_int_m2(0,l2+1,l3,0,m2,m3)+(l2+m2)*J_int_m2(0,l2-1,l3,0,m2,m3)));
      }
      if(l1>=m1+2)
         return sign*((1./(double(l1-m1)*double(l1-m1-1)))*(double(2*l1-1)*gaunt_formula(l1-1,l2,l3,m1+1,m2,m3)
               +double(l1+m1)*(l1+m1-1)*J_int_m2(l1-2,l2,l3,m1,m2,m3)));
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
         return -sign*(2*l2-1)*gaunt_formula(l1,l2-1,l3,m1,m2-1,m3);
      else if(l2 == m2+1 && l2>=2)
         return -sign*(2*m2+1)*gaunt_formula(l1,m2,l3,m1,m2-1,m3);
      else if(l2==1 && m2==0)
      {
         return sign*(1./(2.*l3+1.))*((l3-m3+1)*J_int_m2(0,0,l3+1,0,0,m3)+(l3+m3)*J_int_m2(0,0,l3-1,0,0,m3));
      }

      if(l2>=m2+2)
         return sign*((1./(double(l2-m2)*double(l2-m2-1)))*(double(2*l2-1)*gaunt_formula(l1,l2-1,l3,m1,m2+1,m3)
               +double(l2+m2)*(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3)));

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
         return - sign*(2*l3-1)*gaunt_formula(l1,l2,l3-1,m1,m2,m3-1);
      else if(l3 == m3+1 && l3>=2)
         return - sign*(2*m3+1)*gaunt_formula(l1,l2,m3,m1,m2,m3-1);
      else if(l3==1 && m3==0)
         return 0;

      if(l3>=m3+2)
         return  sign*((1./(double(l3-m3)*double(l3-m3-1)))*(double(2*l3-1)*gaunt_formula(l1,l2,l3-1,m1,m2,m3+1)
               +double(l3+m3)*(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3)));
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
   //CHECKED ON MARCH 2 2020
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
   double sign(1);
   if(m1<0)
   {
      sign*=pow(-1,m1);
      m1=-m1;
   }
   if(m2<0)
   {
      sign*=pow(-1,m2);
      m2=-m2;
   }
   if(m3<0)
   {
      sign*=pow(-1,m3);
      m3=-m3;
   }
   if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;
   if(l1==0 || l1-1 < m1+1)
   {
      return -sign*((1./(2.*l1+1))*gaunt_formula(l1+1,l2,l3,m1+1,m2,m3));
   }
   else
   {
      return sign*((1./(2.*l1+1))*(gaunt_formula(l1-1,l2,l3,m1+1,m2,m3)-gaunt_formula(l1+1,l2,l3,m1+1,m2,m3)));
   }
}
double J_int_p1(int l1,int l2,int l3,int m1,int m2,int m3)
{
   //CHECKED ON MARCH 2 2020
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
   double sign(1);
   if(m1<0)
   {
      sign*=pow(-1,m1);
      m1=-m1;
   }
   if(m2<0)
   {
      sign*=pow(-1,m2);
      m2=-m2;
   }
   if(m3<0)
   {
      sign*=pow(-1,m3);
      m3=-m3;
   }
   if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;
   if(l1==0 || l1-1 < m1)
   {
      return sign*((1./(2.*l1+1.))*(double(l1-m1+1.)*gaunt_formula(l1+1,l2,l3,m1,m2,m3)));
   }
   else
   {
      return sign*((1./(2.*l1+1.))*(double(l1-m1+1.)*gaunt_formula(l1+1,l2,l3,m1,m2,m3)+double(l1+m1)*gaunt_formula(l1-1,l2,l3,m1,m2,m3)));
   }
}
double J_int_D(int l1,int l2,int l3,int m1,int m2,int m3)
{
   double sign(1);
   if(m1<0)
   {
      sign*=pow(-1,m1);
      m1=-m1;
   }
   if(m2<0)
   {
      sign*=pow(-1,m2);
      m2=-m2;
   }
   if(m3<0)
   {
      sign*=pow(-1,m3);
      m3=-m3;
   }
   if( l2 == 0 && l3 == 0  )
      return 0;

   else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
      return 0;

   else if(l2 == 0)
      return sign*((double(l3)/double(l3+m3+1))*(double(l3-m3+1)*J_int_m2(l1,l2,l3+1,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2,m3+1))-double(l3+m3)*J_int_m2(l1,l2,l3-1,m1,m2,m3));
   else if(l3 == 0)
      return sign*((double(l2)/double(l2+m2+1))*(double(l2-m2+1)*J_int_m2(l1,l2+1,l3,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2+1,m3))-double(l2+m2)*J_int_m2(l1,l2-1,l3,m1,m2,m3));
   else
      return sign*((double(l3)/double(l3+m3+1))*(double(l3-m3+1)*J_int_m2(l1,l2,l3+1,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2,m3+1))-double(l3+m3)*J_int_m2(l1,l2,l3-1,m1,m2,m3)
            +(double(l2)/double(l2+m2+1))*(double(l2-m2+1)*J_int_m2(l1,l2+1,l3,m1,m2,m3)-gaunt_formula(l1,l2,l3,m1,m2+1,m3))-double(l2+m2)*J_int_m2(l1,l2-1,l3,m1,m2,m3));
}
double J_int_m1_D(int l1,int l2,int l3,int m1,int m2,int m3)
{
   //CHECKED ON MARCH 10 2020
//   std::cout<<(l1+m1+l2+m2+l3+m3)<<std::endl;
   if((l1+m1+l2+m2+l3+m3)%2==0) //If the integrand is odd
      return 0;
   else
   {
      double sign(1);
      if(m1<0)
      {
         sign*=pow(-1,m1);
         m1=-m1;
      }
      if(m2<0)
      {
         sign*=pow(-1,m2);
         m2=-m2;
      }
      if(m3<0)
      {
         sign*=pow(-1,m3);
         m3=-m3;
      }
      if( l2 == 0 && l3 == 0  )
         return 0;

      else if(l1<0 || l2<0 || l3<0 || m1>l1 || m2>l2 || m3>l3)
         return 0;

      else if(l2 == 0)
         return sign*(-1./(2.*l3+1.))*((l3+1)*(l3+m3)*gaunt_formula(l1,l2,l3-1,m1,m2,m3)-l3*(l3-m3+1)*gaunt_formula(l1,l2,l3+1,m1,m2,m3));
      else if(l3 == 0)
         return sign*(-1./(2.*l2+1.))*((l2+1)*(l2+m2)*gaunt_formula(l1,l2-1,l3,m1,m2,m3)-l2*(l2-m2+1)*gaunt_formula(l1,l2+1,l3,m1,m2,m3));
      else
         return sign*((-1./(2.*l2+1.))*((l2+1)*(l2+m2)*gaunt_formula(l1,l2-1,l3,m1,m2,m3)-l2*(l2-m2+1)*gaunt_formula(l1,l2+1,l3,m1,m2,m3))
               +(-1./(2.*l3+1.))*((l3+1)*(l3+m3)*gaunt_formula(l1,l2,l3-1,m1,m2,m3)-l3*(l3-m3+1)*gaunt_formula(l1,l2,l3+1,m1,m2,m3)));
   }
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
   //CHECKED ON MARCH 10 2020
   if((l1+m1+l2+m2+l3+m3+1)%2 == 0)
      return 0;
   double sign(1);

   if(m1<0)
   {
      sign*=pow(-1,m1);
      m1=-m1;
   }
   if(m2<0)
   {
      sign*=pow(-1,m2);
      m2=-m2;
   }
   if(m3<0)
   {
      sign*=pow(-1,m3);
      m3=-m3;
   }
   if((l2 == 0 && l3 == 0))
      return 0;

   else if(l3 == 0)
   {
      return sign*((double(l2+m2)/double(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3))
      -(double(l2)/double(2*l2+1))*(
            (double(l2-m2+1)/double(2*l2+3))*(double(l2-m2+2)*J_int_m2(l1,l2+2,l3,m1,m2,m3)+(l2+m2+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l2+m2)/(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3))));
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
      return sign*((double(l3+m3)/double(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))
      -(double(l3)/double(2*l3+1))*(
            (double(l3-m3+1)/double(2*l3+3))*(double(l3-m3+2)*J_int_m2(l1,l2,l3+2,m1,m2,m3)+(l3+m3+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l3+m3)/(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))));
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
      return sign*((double(l2+m2)/double(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3))
      -(double(l2)/double(2*l2+1))*(
            (double(l2-m2+1)/double(2*l2+3))*(double(l2-m2+2)*J_int_m2(l1,l2+2,l3,m1,m2,m3)+(l2+m2+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l2+m2)/(2*l2-1))*(double(l2-m2)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l2+m2-1)*J_int_m2(l1,l2-2,l3,m1,m2,m3)))
      + (double(l3+m3)/double(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))
      -(double(l3)/double(2*l3+1))*(
            (double(l3-m3+1)/double(2*l3+3))*(double(l3-m3+2)*J_int_m2(l1,l2,l3+2,m1,m2,m3)+(l3+m3+1)*J_int_m2(l1,l2,l3,m1,m2,m3))
            +(double(l3+m3)/(2*l3-1))*(double(l3-m3)*J_int_m2(l1,l2,l3,m1,m2,m3)+double(l3+m3-1)*J_int_m2(l1,l2,l3-2,m1,m2,m3))));
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

   //Compute the prefactor due to negative degree
   double sign(1);
   if(m1<0)
   {
      m1=-m1;
      sign*=(l1-m1)/(l1+m1);
   }
   if(m2<0)
   {
      m2=-m2;
      sign*=(l2-m2)/(l2+m2);
   }
   if(m3<0)
   {
      m3=-m3;
      sign*=(l3-m3)/(l3+m3);
   }

   //Bypassing algorithm by using fixd header files.
   //
//   if(l1 <= 10 && l2<=10 && l3<=10)
//      return sign*GAUNT_COEFF_VAL[l1*(l1+1)/2+m1][l2*(l2+1)/2+m2][l3*(l3+1)/2+m3];
   //end of bypass. Resuming to the routine

   double temp;
   double factor(0);
   double sum(0);
   double G12(0);
   double G123(0);
   double ALP_integ(0);
   double prefactor(0);


   int m12(m1+m2);
   int m123(m12+m3);

   int fac1[MAX_FACTORIAL_PRIME];
   double temp2(1);
   double val;


   // Checking if any of the polynomials is zero
   
   if(l1 <0 || l2<0 || l3<0)
      return 0;

   else if(m1>l1 || m2>l2 || m3>l3)
      return 0;


   //If none of the polynomials are zero, first compute the factor in front of expansion

      temp=1;
      for (int tt=l1-m1+1;tt!=l1+m1+1;tt++)
      {
         temp*=double(tt);
      }
      val=sqrt(temp);

      temp=1;
      for (int tt=l2-m2+1;tt!=l2+m2+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

      temp=1;
      for (int tt=l3-m3+1;tt!=l3+m3+1;tt++)
      {
         temp*=double(tt);
      }
      val*=sqrt(temp);

      prefactor=val;

   // then, compute the expansion of the Legendre polynomials product
   for(int l12=fabs(l1-l2);l12!=l1+l2+1;l12++)
   {
      if((l1+l2+l12)%2!=0 || m12 > l12)
         continue;
      else
      {

         G12=pow(-1.,m12)*(2.*l12+1.)*WignerSymbols::wigner3j(l1,l2,l12,0,0,0)*WignerSymbols::wigner3j(l1,l2,l12,m1,m2,-m12);
         for(int l123=fabs(l12-l3);l123!=l12+l3+1;l123++)
         {
            if((l12+l3+l123)%2!=0 || m123 > l123)
               continue;
            else
            {
               G123=pow(-1,m123)*(2.*l123+1.)*WignerSymbols::wigner3j(l12,l3,l123,0,0,0)*WignerSymbols::wigner3j(l12,l3,l123,m12,m3,-m123);
               temp=1;
               for (int tt=l123-m123+1;tt!=l123+m123+1;tt++)
               {
                  temp/=double(tt);
               }
               val=sqrt(temp);
               factor=val;

               fact_prime_decomposer(int(double(l123-m123)/2),fac1);
               temp=1;
               for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
               {
                  temp*=pow(double(PRIME[i]),double(fac1[i]));
               }

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
               }


               sum+=G12*G123*factor*ALP_integ;
            }
         }
      }
   }
   return prefactor*sum;
}
/*double Dint(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo)
{
   return sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*acos(-1)))*wigner3j(l1,l2,l3,0,0,0,lnfact_memo)*wigner3j(l1,l2,l3,m1,m2,m3,lnfact_memo);
}*/
