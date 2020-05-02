
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
      sum+=pow(X[i]-av_X,2);
   return sqrt(sum);
}
double covariance(double* X,double* Y,int N)
{

   double av_X(average(X,N));
   double av_Y(average(Y,N));

   double sum(0);

   for(int i=0;i!=N;i++)
      sum+=(X[i]-av_X)*(Y[i]-av_Y);

   return sum;
}
double correlation_coefficient(double* X,double* Y,int N)
{

   double sig_X(sigma(X,N));
   double sig_Y(sigma(Y,N));

   return covariance(X,Y,N)/(sig_X*sig_Y);
}
