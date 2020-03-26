double rYlm (int l,int m,double thet,double phi)
{
//   std::cout<<"Entering rYlm with parameters "<<l<<","<<m<<","<<thet<<","<<phi<<std::endl;
      if(m<0)
      {
         return prefactor_rYlm(l,m)*associated_legendre_nonorm(l,-m,cos(thet))*sin(fabs(m)*phi);
      }
      else if(m>0)
      {
         return prefactor_rYlm(l,m)*associated_legendre_nonorm(l,m,cos(thet))*cos(m*phi);
      }
      else
      {
         return prefactor_rYlm(l,m)*associated_legendre_nonorm(l,m,cos(thet));
      }
}
double prefactor_rYlm(int l,int m)
{
   /*
    * THIS FUNCTION COMPUTES THE SPHERICAL HARMONICS NORMALIZATION PREFACTOR 
    * WE EXPLICITELY INCLUDE THE CONDON SHORTLEY PHASE IN THE ECXPRESSION OF THE SPHERICAL HARMIONICS PREFACTOR
    * SO THAT WE DO NOT TAKE IT  TWICE INTO ACCOUNT WHENEVALUATING THE ASSOCIATED LEGENDRE POLYNOMIALS.
    */
   double temp(1);
   double val(0);

   if(fabs(m) > fabs(l))
   {
      std::cout<<"FATAL ERROR IN SPHERICAL HARMONICS PREFACTOR. M>L:"<<m<<">"<<l<<std::endl;
      return 0;
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

      return val;
   }
   else 
   {
      return pow(-1,m) * prefactor_rYlm(l,fabs(m));
   }
}
/*
 * double associated_legendre(unsigned int l,int m,double x)
 * {
 *  return prefactor_rYlm(l,m)*associated_legendre_nonorm(l,m,x);
 * }
 */
double associated_legendre_nonorm(unsigned int l,int m,double x)
{

   /*
    * COMPUTATION OF THE ASSOCIATED LEGENDRE POLYNOMIAL,!!! ___EXCLUDING THE CONDON-SHORTLEY PHASE___ !!!
    * THE POLYMOMIALS ARE COMUPTED USING RECCURENCE RELATIONS.
    *
    * THE FOLLOWING IDENTITIES ARE USED :
    *
    * P_{L}^{M}(1) = P_{L}^{M}(-1) = 0 FOR M != 0
    *
    * P_{L}^{L}(X) = (2L-1)!! * (1-X**2)^{L/2}
    *
    * P_{L}^{M}(X) = [ (2L-1) / (L-M) ] * X * P_{L-1}^{M}(X) -  [ (L+M-1) / (L-M) ] * P_{L-2}^{M}(X)
    *
    */

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
         if(l == m)
            return dfactorial(2*l-1)*pow(1-x*x,double(l)/2);
         else
            return (double(2*l-1)/double(l-m))*x*associated_legendre_nonorm(l-1,m,x) - double((l+m-1)/double(l-m))*associated_legendre_nonorm(l-2,m,x);
      }
      else 
      {
         return associated_legendre_nonorm(l,-m,x);
      }
   }
}
/*
 * double associated_legendre_der(unsigned int l,int m,double x)
 * {
 *  return associated_legendre_nonorm_der(l,m,x);
 * }
 */ 
double associated_legendre_nonorm_der(unsigned int l,int m,double x)
{

   double val(0);

   if( x == 1 )
      return 0;
   else
   {
      if( m == 0 )
      {
         return legendre_der(l,x);
      }
      else if( m > 0 )
      {
         return (l*x*associated_legendre_nonorm(l,m,x)-(l+m)*associated_legendre_nonorm(l-1,m,x))/sqrt(1-x*x);
      }
      else 
      {
         std::cout<<"Error: Negative ml value in associated legendre der function. exit"<<std::endl;
         exit(EXIT_SUCCESS);
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
