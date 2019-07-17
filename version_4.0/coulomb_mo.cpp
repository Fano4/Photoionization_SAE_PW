/*
 * We implement the Confluent hypergeometric function using 
 * its series definition involving the Pocchamer symbol. 
 *
 * The series is truncated so that the variation of the series
 * value with each term is less than 1e-12
 *
 * M(a,b,z)=SUM_n^inf( A_n )
 * A_n = (a)_n z^n / ( (b)_n * n! ) 
 *     = A_(n-1) * z * ( a + n - 1) / ( n * (b + n - 1) )
 * A_0 = 1
 * ERROR : S_n = A_n - A_(n-1) <= 1e-12
 *
 */
double M_hyperg(double a,double b,double z)
{
   double sum(0);
   double S(1);
   double A(1);
   int n(0);

   while( fabs(S) > 1e-12 )
   {
      sum += A;

      S = -A;
      
      n++;

      A = A * z * ( a + n - 1) / ( n * ( b + n - 1) );
      
      S += A;

   }

   return sum;

}
/*
 * The Boys function is defined recurrently, from the expression
 * the first kind regularized confluent hypergeometric function.
 *
 * We use the downward recurrence as it is more precise.
 *
 * F_(n-1)(z) = ( 2 * z * F_n (z) + exp(-z) ) / ( 2 * n - 1 )
 */
double boys(int n_max,int n,double z, int boys_nmax=nan);
{
   if ( n < n_max)
      return ( 2 * z * boys( n_max , n+1 , z , boys_nmax) + exp(-z) ) / ( 2 * n + 1 );
   else
   {
      if ( isnan(boys_max) )
      {
         boys_nmax = M_hyperg( n_max + 0.5 , n_max +1.5 , -z ) / ( 2 * n_max + 1 );
         return boys_nmax;
      }

      else 
         return boys_nmax;
   }
}
/*
 * Hermite Coulomb function is defined as the integral of the Hermite 
 * Gaussian corresponding to the product between cartesian basis functions
 */
double hermite_coulomb( int n, int t, int u, int v, double p, double x,double y,double z,n_max=nan,boys_nmax=nan)
{
   if(isnan(n_max))
      n_max=t+u+v;

   r=sqrt(x*x+y*y+z*z);

   if( t != 0 )
   {
      if( t > 1 )
         return ( t - 1 ) * hermite_coulomb( n + 1 , t - 2 , u , v , p , x , y , z , n_max , boys_nmax ) + x * hermite_coulomb( n + 1 , t - 1 , u , v , p , x , y , z , n_max , boys_nmax );
      else
         return x * hermite_coulomb( n + 1 , t - 1 , u , v , p , x , y , z , n_max , boys_nmax );

   }

   else if( u != 0 )
   {
      if( u > 1 )
         return ( u - 1 ) * hermite_coulomb( n + 1 , t , u - 2 , v , p , x , y , z , n_max , boys_nmax ) + y * hermite_coulomb( n + 1 , t , u - 1 , v , p , x , y , z , n_max , boys_nmax );
      else
         return y * hermite_coulomb( n + 1 , t , u - 1 , v , p , x , y , z , n_max , boys_nmax );

   }

   else if( v != 0 )
   {
      if( v > 1 )
         return ( v - 1 ) * hermite_coulomb( n + 1 , t , u , v - 2 , p , x , y , z , n_max , boys_nmax ) + z * hermite_coulomb( n + 1 , t , u , v - 2 , p , x , y , z , n_max , boys_nmax );
      else
         return z * hermite_coulomb( n + 1 , t , u , v - 1 , p , x , y , z , n_max , boys_nmax );

   }

   else
   {
      return  pow( -2 * p , n ) * boys( n_max , n , p * r * r , boys_nmax);
   }
}

double E( int i , int j , int t , double xa , double xb , double m , double p )
{
    if( t < 0 || t > i + j)
       return 0;
    else if( t == 0 && i == 0 && j == 0 )
       return exp( -m * pow( xa - xb , 2 ) );
    else if ( t == 0 )
    {
       if( i > 0 )
          return ( xp - xb ) * E( i - 1 , j , 0 , xa , xb , xp , m , p ) + E( i - 1 , j , 1 , xa , xb , xp , m , p ) 
       if( j > 0 )
          return ( xp - xa ) * E( i , j - 1 , 0 , xa , xb , xp , m , p ) + E( i , j - 1 , 1 , xa , xb , xp , m , p ) 
    }
    else
       return ( 1 / ( 2 * p * t ) ) * ( i * E( i - 1 , j , t - 1 , xa , xb , xp , m , p ) + j * E( i , j - 1 , t - 1 , xa , xb , xp , m , p ) )

}

double mono_gauss_prod(int i,int j,int k,int l,int m,int n,double p,double q,double xa,double ya,double za,double xb,double yb,double zb, double xc,double yc,double zc)
{
   double result(0.0);
   double xp(p*xa+q*xb)/(p+q);
   double yp(p*ya+q*yb)/(p+q);
   double zp(p*za+q*zb)/(p+q);

   for(int t=0;t!=i+l+1+1;t++)
   {
      for(int u=0;u!=j+m+1+1;u++)
      {
         for(int v=0;v!=k+n+1+1;v++)
         {
            result+=E(i,l,t,xa,xb,xp,p*q/(p+q),p+q)*E(j,m,u,ya,yb,yp,p*q/(p+q),p+q)*E(k,n,v,za,zb,zp,p*q/(p+q),p+q)*hermite_coulomb(0,t,u,v,p+q, xp-xc,yp-yc,zp-zc)
         }
      }
   }
   result*=2*acos(-1)/(p+q);

   return result;
}

double bi_gauss_prod(int i,int j,int k,int l,int m,int n,int ii,int jj, kk,int ll,int mm,int nn ,double p,double q,double pp,double qq,double xa,double ya,double za,double xb,double yb,double zb)
{
   double result(0.0);
   double xp(p*xa+q*xb)/(p+q);
   double yp(p*ya+q*yb)/(p+q);
   double zp(p*za+q*zb)/(p+q);
   double xq(pp*xa+qq*xb)/(pp+qq);
   double yq(pp*ya+qq*yb)/(pp+qq);
   double zq(pp*za+qq*zb)/(pp+qq);
   double alpha((p+q)*(pp+qq)/((p+q)+(pp+qq)));
   double temp(0);

   for(int t=0;t!=i+l+1+1;t++)
   {
      for(int u=0;u!=j+m+1+1;u++)
      {
         for(int v=0;v!=k+n+1+1;v++)
         {
            temp=0;
            for(int tt=0;tt!=ii+ll+1+1;tt++)
            {
               for(int uu=0;uu!=jj+mm+1+1;uu++)
               {
                  for(int vv=0;vv!=kk+nn+1+1;vv++)
                  {
                     temp+=pow(-1,tt+uu+vv)*E(ii,ll,tt,xa,xb,xq,pp*qq/(pp+qq),pp+qq)*E(jj,mm,uu,ya,yb,yq,pp*qq/(pp+qq),pp+qq)*E(kk,nn,vv,za,zb,zq,pp*qq/(pp+qq),pp+qq)*hermite_coulomb(0,tt+t,uu+u,vv+v,alpha, xp-xq,yp-yq,zp-zq)
                  }
               }
            }
            temp*=E(i,l,t,xa,xb,xp,p*q/(p+q),p+q)*E(j,m,u,ya,yb,yp,p*q/(p+q),p+q)*E(k,n,v,za,zb,zp,p*q/(p+q),p+q)
               result+=temp;
         }
      }
   }
   result*=2*pow(acos(-1),2.5)/((p+q)*(pp+qq)*sqrt((p+q)+(pp+qq)));

   return result;
}
