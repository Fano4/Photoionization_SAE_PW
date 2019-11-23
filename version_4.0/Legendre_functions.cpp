double associated_legendre(unsigned int l,int m,double x,double* lnfact_memo)
{
   int sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));

   if(fabs(m) > fabs(l))
   {
      std::cout<<"FATAL ERROR IN ASSOCIATED LEGENDRE COMPUTATION. M>L:"<<m<<">"<<l<<std::endl;
      return 0;
   }
   else if(fabs(x)==1 && m!=0)
      return 0;
   else
   {
      if(m == 0)
      {
         return sqrt((2*l+1)/(4*acos(-1)))*legendre(l,x);
      }
      else if(m > 0)
      {
         return sign * sqrt((2*l+1) * exp(0.5*(ln_factorial(l-m,lnfact_memo) - ln_factorial(l+m,lnfact_memo))) / (4*acos(-1)))
            * ((l-m+1) * x * associated_legendre_nonorm(l,m-1,x) - (l+m-1) * associated_legendre_nonorm(l-1,m-1,x)) / sqrt(1-x*x);
      }
      else 
      {
         return sign * exp(0.5*(ln_factorial(l+m,lnfact_memo)-ln_factorial(l-m,lnfact_memo)))*associated_legendre(l,-m,x,lnfact_memo);
      }
   }
}
double associated_legendre_nonorm(unsigned int l,int m,double x)
{

   int sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));

   if(fabs(m) > fabs(l))
   {
      std::cout<<"FATAL ERROR IN ASSOCIATED LEGENDRE COMPUTATION. M>L:"<<m<<">"<<l<<std::endl;
      return 0;
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
         return sign * ((l-m+1)*x*associated_legendre_nonorm(l,m-1,x)-(l+m-1)*associated_legendre_nonorm(l-1,m-1,x))/sqrt(1-x*x);
      }
      else 
      {
         return sign * associated_legendre_nonorm(l,-m,x);
      }
   }
}
double associated_legendre_der(unsigned int l,int m,double x,double* lnfact_memo)
{

   int sign(-bool( m % 2 != 0 ) + bool( m % 2 == 0 ));

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
         return -sign * sqrt((2*l+1)*exp(ln_factorial(l-m,lnfact_memo)-ln_factorial(l+m,lnfact_memo))/(4*acos(-1)))
            *(l*x*associated_legendre_nonorm(l,m,x)-(l+m)*associated_legendre_nonorm(l-1,m,x))/sqrt(1-x*x);
      }
      else 
      {
         return sign * exp(ln_factorial(l+m,lnfact_memo)-ln_factorial(l-m,lnfact_memo))*associated_legendre_der(l,-m,x,lnfact_memo);
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
