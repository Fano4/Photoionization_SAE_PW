double associated_legendre(unsigned int l,int m,double x)
{
   double sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));
   double val(1);
   double temp(0);

   if(fabs(m) > fabs(l))
   {
      return 0;
//      std::cout<<"FATAL ERROR IN ASSOCIATED LEGENDRE COMPUTATION. M>L:"<<m<<">"<<l<<std::endl;
//      exit(EXIT_SUCCESS);
   }
   else if(fabs(x)==1 && m!=0)
      return 0;
   else
   {
      if(m == 0)
      {
//          std::cout<<"P_"<<l<<"("<<x<<") : "<<sqrt((2*l+1)/(4*acos(-1)))<<" * "<<legendre(l,x)<<std::endl;
         return sqrt((2*l+1)/(4*acos(-1)))*legendre(l,x);
      }
      else if(m > 0)
      {
         val=sign * sqrt((2*l+1)/ (4*acos(-1)));
         temp=1;
         for (int tt=l-m+1;tt!=l+m+1;tt++)
         {
            temp/=double(tt);
         }
//         std::cout<<"sqrt("<<l<<"!-"<<m<<"!)/("<<l<<"!+"<<m<<"!) = "<<sqrt(temp)<<std::endl;
         val*=sqrt(temp);
         //std::cout<<"===="<<val<<std::endl;
//         std::cout<<" From associated Legendre, calling associated legendre nonorm with parameters "<<l<<","<<m-1<<";"<<l-1<<","<<m-1<<std::endl;
         return val * ((l-m+1) * x * associated_legendre_nonorm(l,m-1,x) - (l+m-1) * associated_legendre_nonorm(l-1,m-1,x)) / sqrt(1-x*x);
      }
      else 
      {
         temp=1;
         for (int tt=l-fabs(m)+1;tt!=l+fabs(m)+1;tt++)
         {
            temp/=tt;
         }
         return sign * temp * associated_legendre(l,-m,x);
      }
   }
}
double associated_legendre_nonorm(unsigned int l,int m,double x)
{

   int sign(1);//(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));

   if(fabs(m) > fabs(l))
   {
      return 0;
//      std::cout<<"FATAL ERROR IN ASSOCIATED LEGENDRE COMPUTATION LINE 52. M>L:"<<m<<">"<<l<<std::endl;
//      exit(EXIT_SUCCESS);
   }
   else if(fabs(x)==1 && m!=0)
      return 0;
   else
   {
      if(m == 0)
      {
         return sign * legendre(l,x);
      }
      else if(m > 0)
      {
//         std::cout<<" Recalling associated legendre nonorm with parameters "<<l<<","<<m-1<<";"<<l-1<<","<<m-1<<std::endl;
         return sign * ((l-m+1)*x*associated_legendre_nonorm(l,m-1,x)-(l+m-1)*associated_legendre_nonorm(l-1,m-1,x))/sqrt(1-x*x);
      }
      else 
      {
//         std::cout<<" From negative m part, Recalling associated legendre nonorm with parameters "<<l<<","<<-m<<std::endl;
         return sign * associated_legendre_nonorm(l,-m,x);
      }
   }
}
double associated_legendre_der(unsigned int l,int m,double x)
{

   int sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));
   double temp(0);
   double val(0);

   if( x == 1 )
      return 0;
   else
   {
      if(m == 0)
      {
         return sqrt((2*l+1)/(4*acos(-1)))*legendre_der(l,x);
      }
      else if(m > 0)
      {
         val=sign * sqrt((2*l+1)/ (4*acos(-1)));
         temp=1;
         for (int tt=l-m+1;tt!=l+m+1;tt++)
         {
            temp/=tt;
         }
         val*=sqrt(temp);
//         std::cout<<" From associated Legendre der, calling associated legendre nonorm with parameters "<<l<<","<<m<<";"<<l-1<<","<<m<<std::endl;
         return -sign * val
//         return -sign * sqrt((2*l+1)*exp(ln_factorial(l-m)-ln_factorial(l+m))/(4*acos(-1)))
            *(l*x*associated_legendre_nonorm(l,m,x)-(l+m)*associated_legendre_nonorm(l-1,m,x))/sqrt(1-x*x);
      }
      else 
      {
         std::cout<<"Legendre function line 95"<<std::endl;
         exit(EXIT_SUCCESS);
//         return sign * exp(ln_factorial(l+m)-ln_factorial(l-m))*associated_legendre_der(l,-m,x);
         //return sign * double(factorial(l+m,fact_memo))*associated_legendre_der(l,-m,x,fact_memo)/double(factorial(l-m,fact_memo));
      }
   }
}
double legendre(unsigned int l,double x)
{
   switch (l)
   {
      case 0:
         return 1;
      case 1:
         return x;
      default:
         return ((2*l-1)*x*legendre(l-1,x)-(l-1)*legendre(l-2,x))/l; 
   }
}
double legendre_der(unsigned int l,double x)
{
   if(x==1)
   {
      return 0;
   }
   else
   {
      switch (l)
      {
         case 0:
            return 0;
         case 1:
            return sqrt(1-x*x);
         default:
            return (l/sqrt(1-x*x))*(x*legendre(l,x)-legendre(l-1,x)); 
      }
   }
}
