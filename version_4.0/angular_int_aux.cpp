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
   double temp;
   double factor(0);
   double sum(0);
   double G12(0);
   double G123(0);
   double ALP_integ(0);
   double prefactor(0);

   // Sorting, sign flips and normalization
   
   Jint_sort_indices(&l1,&l2,&l3,&m1,&m2,&m3); // rearrange to get l1 < l2 < l3
   prefactor=Jint_signflip_renormalize(l1,l2,l3,&m1,&m2,&m3); // Flip the sign of negative m's and renormalize accordingly
   prefactor*=Jint_normalize(l1,l2,l3,m1,m2,m3); // Normalize the product of ALPs

   // Check for zero integrand and for special cases
   if( Jint_special_cases(l1,l2,l3,m1,m2,m3,&temp) ) 
      return prefactor*temp;

   // If the integrand is not special, compute the expansion of the Legendre polynomials product

   //Declare intermediary m's
   
   int m12(m1+m2);
   int m123(m12+m3);

   sum=0;
   for(int l12=abs(l1-l2);l12!=l1+l2+1;l12++)
   {

      if( (l1+l2+l12) % 2 !=0 || m12 > l12) //Here, the selection rules for Wigner 3J symbols apply on each term
         continue;
      else
      {

         G12=pow(-1.,m12)*(2.*l12+1.)*WignerSymbols::wigner3j(l1,l2,l12,0,0,0)*WignerSymbols::wigner3j(l1,l2,l12,m1,m2,-m12);

         for(int l123=fabs(l12-l3);l123!=l12+l3+1;l123++)
         {
            if( (l12+l3+l123) % 2 != 0 || m123 > l123) // Wigner 3J selection rule
               continue;
            else
            {
               G123=pow(-1,m123)*(2.*l123+1.)*WignerSymbols::wigner3j(l12,l3,l123,0,0,0)*WignerSymbols::wigner3j(l12,l3,l123,m12,m3,-m123);

               temp=1;
               for (int tt = l123-m123+1;tt!=l123+m123+1;tt++)
               {
                  temp/=double(tt);
               }
               factor=sqrt(temp);

               sum+=G12*G123*factor*ALP_integral(l123,m123);
            }
         }
      }
   }
   return prefactor*sum;
}
double ALP_integral(int l,int m)
{

   if( (l+m) % 2 != 0)//Check for parity
      return 0;

   else if( l == 0 ) //Check for zero degree
      return 2;
   else
   {
      int fac1[MAX_FACTORIAL_PRIME];
      double temp;
      fact_prime_decomposer(int(double(l-m)/2),fac1);
      temp=1;
      for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
      {
          temp*=pow(double(PRIME[i]),double(fac1[i]));
      }

      return (pow(-1,m)+pow(-1,l))
         *pow(2.,m-2)*m
         *gamma_int_or_half(double(l)/2.)
         *gamma_int_or_half(double(l+m+1)/2.)
         /(temp*gamma_int_or_half(double(l+3)/2.));
   }
}
bool Jint_special_cases(int l1,int l2,int l3,int m1,int m2,int m3,double* result)
{
   int delta(m2+(m1-m2)*bool(m2>=m1));
   double temp;
   double prefactor(0);

   // Checking if any of the polynomials is zero
   
   if(l1 <0 || l2<0 || l3<0) //Negative degree should not occur
   {
      *result=0;
      return 1;
   }
   else if(m1>l1 || m2>l2 || m3>l3) // These ensure the degree of the polynomial is at least zero
   {
      *result=0;
      return 1;
   }

   else if( (l1+l2+l3) % 2 != (m1+m2+m3) % 2 ) // The integrand should be even,meaning that the sums of orders and degree have same parity
   {
      *result=0;
      return 1;
   }

   //If none of the polynomials is zero and the integrand is even,
   //
   // Sorting, sign flips and normalization
   
   Jint_sort_indices(&l1,&l2,&l3,&m1,&m2,&m3); // rearrange to get l1 < l2 < l3
   prefactor=Jint_signflip_renormalize(l1,l2,l3,&m1,&m2,&m3); // Flip the sign of negative m's and renormalize accordingly
   prefactor*=Jint_normalize(l1,l2,l3,m1,m2,m3); // Normalize the product of ALPs
   
   //Then check for special cases

   if( l1 == 0 ) //Overlap integral of two ALPs
   {
      temp=1;
      for (int tt = l2-m2+1;tt!=l2+m2+1;tt++)
      {
          temp*=double(tt);
      }
      *result = bool( l2 == l3) * 2 * temp / (2 * l2 + 1);
      return 1;
   }
   else if(m1+m2==m3)
   {
      *result = 2.*pow(-1,m3)*prefactor*WignerSymbols::wigner3j(l1,l2,l3,0,0,0)*WignerSymbols::wigner3j(l1,l2,l3,m1,m2,-m3);
      return 1;
   }

   else if(abs(m1+m2)==m3)
   {
      *result = 2.*pow(-1,delta-m1+m2)*prefactor*WignerSymbols::wigner3j(l1,l2,l3,0,0,0)*WignerSymbols::wigner3j(l1,l2,l3,-m1,m2,m1-m2);
      return 1;
   }
   else
      return 0;

}
double Jint_normalize(int l1,int l2,int l3,int m1,int m2,int m3)
{
   double temp;
   double val(1);

      temp=1;
      for (int tt=l1-m1+1;tt!=l1+m1+1;tt++)
      {
         temp*=pow(double(tt),1);
      }
      val*=temp;

      temp=1;
      for (int tt=l2-m2+1;tt!=l2+m2+1;tt++)
      {
         temp*=pow(double(tt),1);
      }
      val*=temp;

      temp=1;
      for (int tt=l3-m3+1;tt!=l3+m3+1;tt++)
      {
         temp*=pow(double(tt),1);
      }
      val*=sqrt(temp);

      return val;
}
double Jint_signflip_renormalize(int l1,int l2,int l3,int* m1,int* m2,int* m3)
{

   bool sgnm1( bool ( *m1 < 0 ) );
   bool sgnm2( bool ( *m2 < 0 ) );
   bool sgnm3( bool ( *m3 < 0 ) );
   double temp;
   double val(1);

   //Flip the sign of m's if they are negative
   if(sgnm1)
      *m1=-*m1;
   if(sgnm2)
      *m2=-*m2;
   if(sgnm2)
      *m2=-*m2;

      temp=1;
      if(sgnm1)
      {
          for (int tt=l1-*m1+1;tt!=l1+*m1+1;tt++)
              temp*=pow(double(tt),-1);
      }
      val*=temp;

      temp=1;
      if(sgnm2)
      {
          for (int tt=l2-*m2+1;tt!=l2+*m2+1;tt++)
              temp*=pow(double(tt),-1);
      }
      val*=temp;

      temp=1;
      if(sgnm3)
      {
          for (int tt=l3-*m3+1;tt!=l3+*m3+1;tt++)
              temp*=pow(double(tt),-1);
      }
      val*=temp;

      return pow(-1,sgnm1**m1+sgnm2**m2+sgnm3**m3)*val;

}
void Jint_sort_indices(int* l1,int* l2,int* l3,int* m1,int* m2,int* m3)
{
   //rearrange the array for getting l1<l2<l3
   int tmpl;
   int tmpm;
   if(*l1 > *l2)
   {
      tmpl=*l1;
      tmpm=*m1;
      *l1=*l2;
      *m1=*m2;
      *l2=tmpl;
      *m2=tmpl;
   }
   if(*l1 > *l3)
   {
      tmpl=*l1;
      tmpm=*m1;
      *l1=*l3;
      *m1=*m3;
      *l3=tmpl;
      *m3=tmpl;
   }
   if(*l2 > *l3)
   {
      tmpl=*l2;
      tmpm=*m2;
      *l2=*l3;
      *m2=*m3;
      *l3=tmpl;
      *m3=tmpl;
   }

   return;
}
