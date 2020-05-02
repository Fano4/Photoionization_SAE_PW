
double average(double* X,int N)
{
   double sum(0);

   for(int i=0;i!=N;i++)
   {
      sum+=X[i]/N;
   }
   return sum;
}
double sigma(double* X,int N)
{
   double sum(0);
   double av_X=average(X,N);

   for(int i=0;i!=N;i++)
      sum+=pow(X[i]-av_X,2)/N;
   return sqrt(sum);
}
double covariance(double* X,double* Y,int N)
{

   double av_X(average(X,N));
   double av_Y(average(Y,N));

   double sum(0);

   for(int i=0;i!=N;i++)
      sum+=(X[i]-av_X)*(Y[i]-av_Y)/N;

   return sum;
}
double correlation_coefficient(double* X,double* Y,int N)
{

   double sig_X(sigma(X,N));
   double sig_Y(sigma(Y,N));

   return covariance(X,Y,N)/(sig_X*sig_Y);
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
