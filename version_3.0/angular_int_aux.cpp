
std::complex<double> bessel_cont_angleint(int jl,int jml,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,std::complex<double> **result)
{
   /*
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
    */
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
double azim_integ(int m1,int m2,int m3)
{
   using namespace std;
   double sum(0);

   if(m1 ==0 && m2 ==0 && m3 ==0)
      return 2*acos(-1);
   else if(m1<0 && m2 < 0)
      sum=0.5*(bool(m3==fabs(m1-m2))-bool(m3==fabs(m1+m2)));
   else if(m1<0 && m2>=0)
      sum=0.5*(bool(m3==-fabs(m1-m2))+bool(m3==-fabs(m1+m2)));
   else if(m1>=0 && m2<0)
      sum=0.5*(bool(m3==-fabs(m1-m2))-bool(m3==-fabs(m1+m2)));
   else
      sum=0.5*(bool(m3==fabs(m1-m2))+bool(m3==fabs(m1+m2)));

//   std::cout<<"azimint probe : "<<sum<<std::endl;
//   std::cout<<m1<<","<<m2<<","<<m3<<std::endl;

   return 2*acos(-1)*sum;
}
double rYlm (int l,int m,double thet,double phi)
{
//   std::cout<<"Entering rYlm with parameters "<<l<<","<<m<<","<<thet<<","<<phi<<std::endl;
   if(m<0)
   {
      return sqrt(2)*associated_legendre(l,-m,cos(thet))*sin(-m*phi);
   }
   else if(m>0)
   {
      return sqrt(2)*associated_legendre(l,m,cos(thet))*cos(m*phi);
   }
   else
   {
      return associated_legendre(l,m,cos(thet));
   }
}
double prefactor_rYlm(int l,int m)
{
   int sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));

   if(fabs(m) > fabs(l))
   {
      std::cout<<"FATAL ERROR IN SPHERICAL HARMONICS PREFACTOR. M>L:"<<m<<">"<<l<<std::endl;
      exit(EXIT_SUCCESS);
   }

   if(m == 0)
   {
      return sqrt((2*l+1)/(4*acos(-1)));
   }
   else if(m > 0)
   {
      return sign * sqrt((2*l+1) * factorial(l-m) / (4*acos(-1) * factorial(l+m)));
   }
   else 
   {
      return sign * sqrt((2*l+1) * factorial(l+m) / (4*acos(-1) * factorial(l-m)));
//      return sign * factorial(l+m)*prefactor_rYlm(l,-m)/factorial(l-m);
   }
}
double gaunt_formula(int l1,int l2,int l3,int m1,int m2,int m3)
{
   //DETERMINE THE LARGEST ORDER
   int temp(0);
   double dtemp(0);
   int u(0);
   int v(0);
   int w(0);
   int l(0);
   int m(0);
   int n(0);
   if(m2>m1)
   {
      if(m3>m2) // m3 > m2 > m1
      {
         u=m3;v=m2;w=m1;
         l=l3;m=l2;n=l1;
      }
      else if (m3 >m1) // m2 > m3 > m1 
      {
         u=m2;v=m3;w=m1;
         l=l2;m=l3;n=l1;
      }
   }
   else if(m3 > m1) //m3 > m1 > m2
   {
         u=m3;v=m1;w=m2;
         l=l3;m=l1;n=l2;
   }
   else
   {
      if(m3 > m2) // m1 > m3 > m2
      {
         u=m1;v=m3;w=m2;
         l=l1;m=l3;n=l2;
      }
      else // m1 > m2 > m3
      {
         u=m1;v=m2;w=m3;
         l=l1;m=l2;n=l3;
      }
   }
   //CHECK THE ORDER OF THE DEGREES
   if(m<n)
   {
      temp=m;
      m=n;
      n=temp;
      temp=v;
      v=w;
      w=temp;
   }

   if(fabs(m1) > l1 || fabs(m2) > l2 || fabs(m3) > l3)
   {
//      std::cout<<"FATAL ERROR: ORDER LARGER THAN DEGREE IN GAUNT FORMULA"<<std::endl;
//      std::cout<<m1<<"<"<<l1<<","<<m2<<"<"<<l2<<","<<m3<<"<"<<l3<<std::endl;
      return 0;
   }
   //CHECK THAT THE PARAMETERS RESPECT THE CONDITIONS
   else if( m1 < 0 || m2 < 0 || m3 < 0)
   {
      std::cout<<"WARNING !! NEGATIVE ORDERS IN GAUNT FORMULA !! "<<std::endl;
      return 0;
   }
   else if (l1 < 0 || l2 < 0 || l3 < 0)
      return 0;
   else if( (l+m+n) % 2 != 0)
      return 0;
   else if( l3 > l1 + l2 || l3 < fabs(l1- l2) )
      return 0;
   else if (u != v + w)
      return 0;
   else
   {
   //IF THE INTEGRAL IS NOT DEFINED ZERO, COMPUTE THE ACTUAL VALUE
//   std::cout<<"gaunt probe...";
   int s((l+m+n)/2);
   int p(std::max(0,n-m-u));
   int q(std::min(m+n-u,std::min(l-u,n-w)));

   double sum(0);

//   std::cout<<u<<","<<v<<","<<w<<","<<l<<","<<m<<","<<n<<std::endl<<s<<","<<p<<","<<q<<std::endl;
   for(int t=p;t<q+1;t++)
   {
//      std::cout<<t<<std::endl;
      sum+=pow(-1,t)*factorial(l+u+t)*factorial(m+n-u-t)/(factorial(t)*factorial(l-u-t)*factorial(m-n+u+t)*factorial(n-w-t));
   }
//   std::cout<<m-v<<","<<s-l<<","<<s-m<<","<<s-n<<","<<2*s+1<<std::endl;
   dtemp=2*pow(-1,s-m-w)*sum*factorial(m+v)*factorial(n+w)*factorial(2*s-2*n)*factorial(s)/( factorial(m-v)*factorial(s-l)*factorial(s-m)*factorial(s-n)*factorial(2*s+1));
//   std::cout<<"Gaunt val"<<dtemp<<std::endl;
   if(isnan(dtemp))
       std::cout<<" ERROR ! GAUNT FUNCTION IS NAN"<<std::endl;
   return dtemp;
   }
}
double Dint(int l1,int l2,int l3,int m1,int m2,int m3)
{
   return sqrt(((2*l1+1)*(2*l2+1)*(2*l3+1))/(4*acos(-1)))*wigner3j(l1,l2,l3,0,0,0)*wigner3j(l1,l2,l3,m1,m2,m3);
}
