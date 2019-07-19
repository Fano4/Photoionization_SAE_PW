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
double boys(int n_max,int n,double z, double boys_nmax)
{
   if ( n < n_max)
      return ( 2 * z * boys( n_max , n+1 , z , boys_nmax) + exp(-z) ) / ( 2 * n + 1 );
   else
   {
      if ( isnan(boys_nmax) )
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
double hermite_coulomb( int n, int t, int u, int v, double p, double x,double y,double z,int n_max,double boys_nmax)
{

   if(n_max==-1)
      n_max=t+u+v;

   double r(sqrt(x*x+y*y+z*z));

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

double E( int i , int j , int t , double xa , double xb,double xp , double m , double p )
{
    if( t < 0 || t > i + j)
       return 0;
    else if( t == 0 && i == 0 && j == 0 )
       return exp( -m * pow( xa - xb , 2 ) );
    else if ( t == 0 )
    {
       if( i > 0 )
          return ( xp - xb ) * E( i - 1 , j , 0 , xa , xb , xp , m , p ) + E( i - 1 , j , 1 , xa , xb , xp , m , p ) ;
       if( j > 0 )
          return ( xp - xa ) * E( i , j - 1 , 0 , xa , xb , xp , m , p ) + E( i , j - 1 , 1 , xa , xb , xp , m , p ) ;
    }
    else
       return ( 1 / ( 2 * p * t ) ) * ( i * E( i - 1 , j , t - 1 , xa , xb , xp , m , p ) + j * E( i , j - 1 , t - 1 , xa , xb , xp , m , p ) );

}

double mono_gauss_prod(int lx1,int ly1,int lz1,int lx2,int ly2,int lz2,double d1,double d2,double x1,double y1,double z1,double x2,double y2,double z2, double xc,double yc,double zc)
{
   double result(0.0);
   double p(d1+d2);
   double mu((d1*d2)/(d1+d2));
   double xp((d1*x1+d2*x2)/(d1+d2));
   double yp((d1*y1+d2*y2)/(d1+d2));
   double zp((d1*z1+d2*z2)/(d1+d2));

   for(int t=0;t!=lx1+lx2+1+1;t++)
   {
      for(int u=0;u!=ly1+ly2+1+1;u++)
      {
         for(int v=0;v!=lz1+lz2+1+1;v++)
         {
            result+=E(lx1,lx2,t,x1,x2,xp,mu,p)*E(ly1,ly2,u,y1,y2,yp,mu,p)*E(lz1,lz2,v,z1,z2,zp,mu,p)*hermite_coulomb(0,t,u,v,p,xp-xc,yp-yc,zp-zc);
         }
      }
   }
   result*=2*acos(-1)/p;

   return result;
}

double bi_gauss_prod(int lxa1,int lya1,int lza1,int lxb1,int lyb1,int lzb1,int lxa2,int lya2,int lza2,int lxb2,int lyb2,int lzb2 ,double da1,double db1,double da2,double db2,double xa1,double ya1,double za1,double xb1,double yb1,double zb1,double xa2,double ya2,double za2,double xb2,double yb2,double zb2)
{
   double result(0.0);
   double p1(da1+db1);
   double p2(da2+db2);
   double mu1(da1*db1/p1);
   double mu2(da2*db2/p2);
   double xp1((da1*xa1+db1*xb1)/(da1+db1));
   double yp1((da1*ya1+db1*yb1)/(da1+db1));
   double zp1((da1*za1+db1*zb1)/(da1+db1));
   double xp2((da2*xa2+db2*xb2)/(da2+db2));
   double yp2((da2*ya2+db2*yb2)/(da2+db2));
   double zp2((da2*za2+db2*zb2)/(da2+db2));
   double alpha(p1*p2/(p1+p2));
   double temp(0);

   for(int t=0;t!=lxa1+lxb1+1+1;t++)
   {
      for(int u=0;u!=lya1+lyb1+1+1;u++)
      {
         for(int v=0;v!=lza1+lzb1+1+1;v++)
         {
            temp=0;
            for(int tt=0;tt!=lxa2+lxb2+1+1;tt++)
            {
               for(int uu=0;uu!=lya2+lyb2+1+1;uu++)
               {
                  for(int vv=0;vv!=lza2+lzb2+1+1;vv++)
                  {
                     temp+=pow(-1,tt+uu+vv)*E(lxa2,lxb2,tt,xa2,xb2,xp2,mu2,p2)*E(lya2,lyb2,uu,ya2,yb2,yp2,mu2,p2)*E(lza2,lzb2,vv,za2,zb2,zp2,mu2,p2)*hermite_coulomb(0,tt+t,uu+u,vv+v,alpha, xp1-xp2,yp1-yp2,zp1-zp2);
                  }
               }
            }
            temp*=E(lxa1,lxb1,t,xa1,xb1,xp1,mu1,p1)*E(lya1,lyb1,u,ya1,yb1,yp1,mu1,p1)*E(lza1,lzb1,v,za1,zb1,zp1,mu1,p1);
            result+=temp;
         }
      }
   }
   result*=2*pow(acos(-1),2.5)/(p1*p2*sqrt(p1+p2));

   return result;
}

double spher_harmo_mono_gauss_int(int l, int lp, int m, int mp,double d,double dp,double xa,double ya,double za,double xb,double yb,double zb,double xc,double yc,double zc)
{
   double result=0;

   int i,j,k,ll,mm,n;
   double C,Cp;
   double vm=0.5*bool(m<0);
   double vmp=0.5*bool(mp<0);

   double N((1/(pow(2,fabs(m))*factorial(l)))*sqrt(2*factorial(l+fabs(m))*factorial(l-fabs(m))/2*bool(m==0)));
   double Np((1/(pow(2,fabs(mp))*factorial(lp)))*sqrt(2*factorial(lp+fabs(mp))*factorial(lp-fabs(mp))/2*bool(mp==0)));

   double Renorm(sqrt(0.5*intplushalf_gamma(1+l)/(pow(2*d,1.5+l)))/pow(acos(-1)/(2*d),0.75));
   double Renormp(sqrt(0.5*intplushalf_gamma(1+lp)/(pow(2*dp,1.5+lp)))/pow(acos(-1)/(2*dp),0.75));
   
   for(int t=0;t!=int((l-fabs(m))/2)+1;t++)
   {
      for(int tp=0;tp!=int((lp-fabs(mp))/2)+1;tp++)
      {
         for(int u=0;u!=t+1;u++)
         {
            for(int up=0;up!=tp+1;up++)
            {
               for(int v=0;v!=int(fabs(m/2)-vm);v++)
               {
                  for(int vp=0;vp!=int(fabs(mp/2)-vmp);vp++)
                  {
                     i=2*t+fabs(m)-2*(u+v);
                     j=2*(u+v);
                     k=l-2*fabs(m);
                     ll=2*tp+fabs(mp)-2*(up+vp);
                     mm=2*(up+vp);
                     n=lp-2*fabs(mp);
                     C=(pow(-1.0,float(t)+float(v)-float(vm)))*(pow(0.25,t))*binomial(l,t)*binomial(l-t,abs(m)+t)*binomial(t,u)*binomial(fabs(m),2*v);
                     Cp=(pow(-1.0,float(tp)+float(vp)-float(vmp)))*(pow(0.25,tp))*binomial(lp,tp)*binomial(lp-tp,abs(mp)+tp)*binomial(tp,up)*binomial(fabs(mp),2*vp);
                     result+=C*Cp*mono_gauss_prod(i,j,k,ll,mm,n,d,dp,xa,ya,za,xb,yb,zb,xc,yc,zc);

                  }
               }
            }
         }
      }
   }
   return result*N*Np*Renorm*Renormp;
}
double spher_harmo_biel_gauss_int(int la1,int lb1, int la2,int lb2, int ma1,int mb1, int ma2,int mb2,double da1,double db1,double da2,double db2,double xa1,double ya1,double za1,double xb1,double yb1,double zb1,double xa2,double ya2,double za2,double xb2,double yb2,double zb2)
{
   double result=0;

   double vma1=0.5*bool(ma1<0);
   double vmb1=0.5*bool(mb1<0);
   double vma2=0.5*bool(ma2<0);
   double vmb2=0.5*bool(mb2<0);

   double Na1=(1/(pow(2,fabs(ma1))*factorial(la1)))*sqrt(2*factorial(la1+fabs(ma1))*factorial(la1-fabs(ma1))/2*bool(ma1==0));
   double Nb1=(1/(pow(2,fabs(mb1))*factorial(lb1)))*sqrt(2*factorial(lb1+fabs(mb1))*factorial(lb1-fabs(mb1))/2*bool(mb1==0));
   double Na2=(1/(pow(2,fabs(ma2))*factorial(la2)))*sqrt(2*factorial(la2+fabs(ma2))*factorial(la2-fabs(ma2))/2*bool(ma2==0));
   double Nb2=(1/(pow(2,fabs(mb2))*factorial(lb2)))*sqrt(2*factorial(lb2+fabs(mb2))*factorial(lb2-fabs(mb2))/2*bool(mb2==0));

   double Renorma1(sqrt(0.5*intplushalf_gamma(1+la1)/(pow(2*da1,1.5+la1)))/pow(acos(-1)/(2*da1),0.75));
   double Renormb1(sqrt(0.5*intplushalf_gamma(1+lb1)/(pow(2*db1,1.5+lb1)))/pow(acos(-1)/(2*db1),0.75));
   double Renorma2(sqrt(0.5*intplushalf_gamma(1+la2)/(pow(2*da2,1.5+la2)))/pow(acos(-1)/(2*da2),0.75));
   double Renormb2(sqrt(0.5*intplushalf_gamma(1+lb2)/(pow(2*db2,1.5+lb2)))/pow(acos(-1)/(2*db2),0.75));

   int lxa1,lxa2,lya1,lya2,lza1,lza2,lxb1,lxb2,lyb1,lyb2,lzb1,lzb2;

   double Ca1,Cb1,Ca2,Cb2;
   
   for(int ta1=0;ta1!=int((la1-fabs(ma1))/2)+1;ta1++)
   {
      for(int tb1=0;tb1!=int((lb1-fabs(mb1))/2)+1;tb1++)
      {
         for(int ta2=0;ta2!=int((la2-fabs(ma2))/2)+1;ta2++)
         {
            for(int tb2=0;tb2!=int((lb2-fabs(mb2))/2)+1;tb2++)
            {
               for(int ua1=0;ua1!=ta1+1;ua1++)
               {
                  for(int ub1=0;ub1!=tb1+1;ub1++)
                  {
                     for(int ua2=0;ua2!=ta2+1;ua2++)
                     {
                        for(int ub2=0;ub2!=tb2+1;ub2++)
                        {
                           for(int va1=0;va1!=int(fabs(ma1/2)-vma1);va1++)
                           {
                              for(int vb1=0;vb1!=int(fabs(mb1/2)-vmb1);vb1++)
                              {
                                 for(int va2=0;va2!=int(fabs(ma2/2)-vma2);va2++)
                                 {
                                    for(int vb2=0;vb2!=int(fabs(mb2/2)-vmb2);vb2++)
                                    {
                                       lxa1=2*ta1+fabs(ma1)-2*(ua1+va1);
                                       lxb1=2*tb1+fabs(mb1)-2*(ub1+vb1);
                                       lxa2=2*ta2+fabs(ma2)-2*(ua2+va2);
                                       lxb2=2*tb2+fabs(mb2)-2*(ub2+vb2);
                                       lya1=2*(ua1+va1);
                                       lyb1=2*(ub1+vb1);
                                       lya2=2*(ua2+va2);
                                       lyb2=2*(ub2+vb2);
                                       lza1=la1-2*fabs(ma1);
                                       lzb1=lb1-2*fabs(mb1);
                                       lza2=la2-2*fabs(ma2);
                                       lzb2=lb2-2*fabs(mb2);
                                       Ca1=(pow(-1.0,float(ta1)+float(va1)-float(vma1)))*(pow(0.25,ta1))*binomial(la1,ta1)*binomial(la1-ta1,abs(ma1)+ta1)*binomial(ta1,ua1)*binomial(fabs(ma1),2*va1);
                                       Cb1=(pow(-1.0,float(tb1)+float(vb1)-float(vmb1)))*(pow(0.25,tb1))*binomial(lb1,tb1)*binomial(lb1-tb1,abs(mb1)+tb1)*binomial(tb1,ub1)*binomial(fabs(mb1),2*vb1);
                                       Ca2=(pow(-1.0,float(ta2)+float(va2)-float(vma2)))*(pow(0.25,ta2))*binomial(la2,ta2)*binomial(la2-ta2,abs(ma2)+ta2)*binomial(ta2,ua2)*binomial(fabs(ma2),2*va2);
                                       Cb2=(pow(-1.0,float(tb2)+float(vb2)-float(vmb2)))*(pow(0.25,tb2))*binomial(lb2,tb2)*binomial(lb2-tb2,abs(mb2)+tb2)*binomial(tb2,ub2)*binomial(fabs(mb2),2*vb2);
                                       result+=Ca1*Cb1*Ca2*Cb2*bi_gauss_prod(lxa1,lya1,lza1,lxb1,lyb1,lzb1,lxa2,lya2,lza2,lxb2,lyb2,lzb2,da1,db1,da2,db2,xa1,ya1,za1,xb1,yb1,zb1,xa2,ya2,za2,xb2,yb2,zb2);

                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return result*Na1*Nb1*Na2*Nb2*Renorma1*Renormb1*Renorma2*Renormb2;
}

bool AO_mono_coulomb(double **AO_mono_coul_mat,int basis_size,int cont_num,int num_of_nucl, double **contraction_coeff_array,double **contraction_zeta_array,double **nucl_spher_pos,int *nucl_basis_func,int **angular_mom_numbers)
{
   double xa(0);
   double ya(0);
   double za(0);
   double xb(0);
   double yb(0);
   double zb(0);
   double xc(0);
   double yc(0);
   double zc(0);

   for(int n=0;n!=num_of_nucl;n++)
   {
      xc=nucl_spher_pos[n][0]*sin(nucl_spher_pos[n][1])*cos(nucl_spher_pos[n][2]);
      yc=nucl_spher_pos[n][0]*sin(nucl_spher_pos[n][1])*sin(nucl_spher_pos[n][2]);
      zc=nucl_spher_pos[n][0]*cos(nucl_spher_pos[n][1]);

      for(int i=0;i!=basis_size;i++)
      {
          xa=nucl_spher_pos[nucl_basis_func[i]][0]*sin(nucl_spher_pos[nucl_basis_func[i]][1])*cos(nucl_spher_pos[nucl_basis_func[i]][2]);
          ya=nucl_spher_pos[nucl_basis_func[i]][0]*sin(nucl_spher_pos[nucl_basis_func[i]][1])*sin(nucl_spher_pos[nucl_basis_func[i]][2]);
          za=nucl_spher_pos[nucl_basis_func[i]][0]*cos(nucl_spher_pos[nucl_basis_func[i]][1]);
         for(int j=0;j!=basis_size;j++)
         {
             xb=nucl_spher_pos[nucl_basis_func[j]][0]*sin(nucl_spher_pos[nucl_basis_func[j]][1])*cos(nucl_spher_pos[nucl_basis_func[j]][2]);
             yb=nucl_spher_pos[nucl_basis_func[j]][0]*sin(nucl_spher_pos[nucl_basis_func[j]][1])*sin(nucl_spher_pos[nucl_basis_func[j]][2]);
             zb=nucl_spher_pos[nucl_basis_func[j]][0]*cos(nucl_spher_pos[nucl_basis_func[j]][1]);
            AO_mono_coul_mat[n][i*basis_size+j]=0;

            for(int p=0;p!=cont_num;p++)
            {
               for(int q=0;q!=cont_num;q++)
               {
                  AO_mono_coul_mat[n][i*basis_size+j]+=contraction_coeff_array[i][p]*contraction_coeff_array[j][q]*spher_harmo_mono_gauss_int(angular_mom_numbers[i][0],angular_mom_numbers[i][1],angular_mom_numbers[j][0],angular_mom_numbers[j][1],contraction_zeta_array[i][p],contraction_zeta_array[j][q],xa, ya, za, xb, yb, zb, xc, yc, zc);
               }
            }
         }
      }
   }
   return 1;
}
bool AO_biel_coulomb(double **AO_biel_coul_mat,int basis_size,int cont_num,int num_of_nucl, double **contraction_coeff_array,double **contraction_zeta_array,double **nucl_spher_pos,int *nucl_basis_func,int **angular_mom_numbers)
{
   double xa1(0);
   double ya1(0);
   double za1(0);
   double xb1(0);
   double yb1(0);
   double zb1(0);
   double xa2(0);
   double ya2(0);
   double za2(0);
   double xb2(0);
   double yb2(0);
   double zb2(0);

   for(int a1=0;a1!=basis_size;a1++)
   {
      xa1=nucl_spher_pos[nucl_basis_func[a1]][0]*sin(nucl_spher_pos[nucl_basis_func[a1]][1])*cos(nucl_spher_pos[nucl_basis_func[a1]][2]);
      ya1=nucl_spher_pos[nucl_basis_func[a1]][0]*sin(nucl_spher_pos[nucl_basis_func[a1]][1])*sin(nucl_spher_pos[nucl_basis_func[a1]][2]);
      za1=nucl_spher_pos[nucl_basis_func[a1]][0]*cos(nucl_spher_pos[nucl_basis_func[a1]][1]);
      for(int b1=0;b1!=basis_size;b1++)
      {
         xb1=nucl_spher_pos[nucl_basis_func[b1]][0]*sin(nucl_spher_pos[nucl_basis_func[b1]][1])*cos(nucl_spher_pos[nucl_basis_func[b1]][2]);
         yb1=nucl_spher_pos[nucl_basis_func[b1]][0]*sin(nucl_spher_pos[nucl_basis_func[b1]][1])*sin(nucl_spher_pos[nucl_basis_func[b1]][2]);
         zb1=nucl_spher_pos[nucl_basis_func[b1]][0]*cos(nucl_spher_pos[nucl_basis_func[b1]][1]);
         for(int a2=0;a2!=basis_size;a2++)
         {
            xa2=nucl_spher_pos[nucl_basis_func[a2]][0]*sin(nucl_spher_pos[nucl_basis_func[a2]][1])*cos(nucl_spher_pos[nucl_basis_func[a2]][2]);
            ya2=nucl_spher_pos[nucl_basis_func[a2]][0]*sin(nucl_spher_pos[nucl_basis_func[a2]][1])*sin(nucl_spher_pos[nucl_basis_func[a2]][2]);
            za2=nucl_spher_pos[nucl_basis_func[a2]][0]*cos(nucl_spher_pos[nucl_basis_func[a2]][1]);
            for(int b2=0;b2!=basis_size;b2++)
            {
               xb2=nucl_spher_pos[nucl_basis_func[b2]][0]*sin(nucl_spher_pos[nucl_basis_func[b2]][1])*cos(nucl_spher_pos[nucl_basis_func[b2]][2]);
               yb2=nucl_spher_pos[nucl_basis_func[b2]][0]*sin(nucl_spher_pos[nucl_basis_func[b2]][1])*sin(nucl_spher_pos[nucl_basis_func[b2]][2]);
               zb2=nucl_spher_pos[nucl_basis_func[b2]][0]*cos(nucl_spher_pos[nucl_basis_func[b2]][1]);

               AO_biel_coul_mat[a1*basis_size+b1][a2*basis_size+b2]=0;

               for(int pa1=0;pa1!=cont_num;pa1++)
               {
                  for(int pb1=0;pb1!=cont_num;pb1++)
                  {
                     for(int pa2=0;pa2!=cont_num;pa2++)
                     {
                        for(int pb2=0;pb2!=cont_num;pb2++)
                        {
                           AO_biel_coul_mat[a1*basis_size+b1][a2*basis_size+b2]+=contraction_coeff_array[a1][pa1]*contraction_coeff_array[b1][pb1]*contraction_coeff_array[a2][pa2]*contraction_coeff_array[b2][pb2]*spher_harmo_biel_gauss_int(angular_mom_numbers[a1][0],angular_mom_numbers[b1][0],angular_mom_numbers[a2][0],angular_mom_numbers[b2][0],angular_mom_numbers[a1][1],angular_mom_numbers[b1][1],angular_mom_numbers[a2][1],angular_mom_numbers[b2][1],contraction_zeta_array[a1][pa1],contraction_zeta_array[b1][pb1],contraction_zeta_array[a2][pa2],contraction_zeta_array[b2][pb2],xa1, ya1, za1, xb1, yb1, zb1,xa2,ya2,za2,xb2,yb2,zb2);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return 1;
}
bool MO_mono_coulomb(double **MO_mono_coul_mat,int n_occ,int n_closed,int basis_size,int cont_num,int num_of_nucl,double *lcao_coeff_array, double **contraction_coeff_array,double **contraction_zeta_array,double **nucl_spher_pos,int *nucl_basis_func,int **angular_mom_numbers)
{
   double **AO_mono_coul_mat=new double*[num_of_nucl];
   double *temp=new double[basis_size*(n_occ+n_closed)];

   for(int n=0;n!=num_of_nucl;n++)
   {
      AO_mono_coul_mat[n]=new double[basis_size*basis_size];
   }

   AO_mono_coulomb(AO_mono_coul_mat,basis_size,cont_num,num_of_nucl,contraction_coeff_array,contraction_zeta_array,nucl_spher_pos,nucl_basis_func,angular_mom_numbers);

   for(int n=0;n!=num_of_nucl;n++)
   {
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,basis_size,n_occ+n_closed,basis_size,1,AO_mono_coul_mat[n],basis_size,lcao_coeff_array,(n_occ+n_closed),0,temp,(n_occ+n_closed));
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n_occ+n_closed,n_occ+n_closed,basis_size,1,lcao_coeff_array,basis_size,temp,n_occ+n_closed,0,MO_mono_coul_mat[n],(n_occ+n_closed));
      /*
      for(int i=0;i!=n_occ+n_closed;i++)
      {
         for(int j=0;j!=n_occ+n_closed;j++)
         {
            MO_mono_coul_mat[n][i*(n_occ+n_closed)+j]=0.0;

            for(int k=0;k!=basis_size;k++)
            {
               for(int l=0;l!=basis_size;l++)
               {
                  MO_mono_coul_mat[n][i*(n_occ+n_closed)+j]+=AO_mono_coul_mat[n][k*basis_size+l]*lcao_coeff_array[i*basis_size+k]*lcao_coeff_array[j*basis_size+l];
               }
            }
         }
      }
      */
   }

   delete [] AO_mono_coul_mat;
   delete [] temp;

   return 1;
}

double MO_biel_coulomb(double **MO_biel_coul_mat,int n_occ,int n_closed,int basis_size,int cont_num,int num_of_nucl,double *lcao_coeff_array, double **contraction_coeff_array,double **contraction_zeta_array,double **nucl_spher_pos,int *nucl_basis_func,int **angular_mom_numbers)
{
   double **AO_biel_coul_mat=new double*[basis_size*basis_size];
   for(int n=0;n!=basis_size*basis_size;n++)
   {
      AO_biel_coul_mat[n]=new double[basis_size*basis_size];
   }

   AO_biel_coulomb(AO_biel_coul_mat,basis_size,cont_num,num_of_nucl,contraction_coeff_array,contraction_zeta_array,nucl_spher_pos,nucl_basis_func,angular_mom_numbers);

      for(int moa1=0;moa1!=n_occ+n_closed;moa1++)
      {
         for(int mob1=0;mob1!=n_occ+n_closed;mob1++)
         {
            for(int moa2=0;moa2!=n_occ+n_closed;moa2++)
            {
               for(int mob2=0;mob2!=n_occ+n_closed;mob2++)
               {
                  MO_biel_coul_mat[moa1*(n_occ+n_closed)+mob1][moa2*(n_occ+n_closed)+mob2]=0.0;

                  for(int aoa1=0;aoa1!=basis_size;aoa1++)
                  {
                     for(int aob1=0;aob1!=basis_size;aob1++)
                     {
                        for(int aoa2=0;aoa2!=basis_size;aoa2++)
                        {
                           for(int aob2=0;aob2!=basis_size;aob2++)
                           {
                              MO_biel_coul_mat[moa1*(n_occ+n_closed)+mob1][moa2*(n_occ+n_closed)+mob2]+=AO_biel_coul_mat[aoa1*basis_size+aob1][aoa2*basis_size+aob2]*lcao_coeff_array[moa1*basis_size+aoa1]*lcao_coeff_array[mob1*basis_size+aob1]*lcao_coeff_array[moa2*basis_size+aoa2]*lcao_coeff_array[mob2*basis_size+aob2];
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   delete [] AO_biel_coul_mat;
   return 1;
}
