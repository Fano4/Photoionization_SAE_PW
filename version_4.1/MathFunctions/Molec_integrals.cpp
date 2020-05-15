#include <cstdlib>
#include <string>
#include <cmath>
#include "mathfunctions.h"

/*bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff)
{
   bool test(0);
   bool test2(0);
   int p;
   int q;

   double *temp=new double[(n_elec_neut)*(n_elec_neut)];

    for (int i=0; i!=n_states_neut; i++)//ELECTRONIC STATE N
    {
        for (int j=0; j!=n_states_cat; j++)//ELECTRONIC STATE K
        {
            for (int k=0; k!=n_occ; k++)//MOLECULAR ORBITAL k COEFF. FOR THE DYSON
            {
                Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]=0;
                test=0;
                
                for (int n=0; n!=ci_size_neut; n++)//   over configurations of the neutral
                {
                    for (int l=0; l!=ci_size_cat; l++)//  over configuration of the cation
                    {
                        test2=0;
                        for(int m=0; m!=n_elec_neut; m++)//Over the electrons of the neutral
                        {
                           if (ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+m]==k && !test2)
                           {
                              test=1;
                              test2=1;
                              continue;
                           }
                           for (int o=0; o!=n_elec_neut-1; o++)//Over the electrons of the cation
                           {
                              p=ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+m];
                              q=ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+o];
                              temp[(n_elec_neut-1)*(m-test2)+o]=overlap[n_occ*p+q]*kronecker_delta(ci_vec_neut[1][n_elec_neut*n+(m-test2)], ci_vec_cat[1][(n_elec_neut-1)*l+o]);;
                           }
                        }
                        if(test2)
                        {
                           Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]+=(ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]*ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]*determinant(temp,(n_elec_neut-1)));
                        }
                    }
                }
                if(!test)
                {
                    Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]=0;
                    std::cout<<std::endl<<"states "<<i<<"  and  "<<j<<"    "<<Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]<<"   MO  "<<k<<std::endl<<"====================================="<<std::endl<<std::endl;
                }
                else
                {
                    std::cout<<std::endl<<"states "<<i<<"  and  "<<j<<"    "<<Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]<<"   MO  "<<k<<std::endl<<"====================================="<<std::endl<<std::endl;
                }
                    
            }
        }
    }


    delete [] temp;

    return 1;
}*/
double bessel_gaussian_poly_integral(int l,double p,double r)
{
   return sqrt(acos(-1))*pow(r,l)*exp(-r*r/(4*p))/(pow(2,2+l)*pow(p,1.5+l));
}
double prim_radial_ovlp(int la,int lb,double zet_a,double zet_b,double r)
{
   return bessel_gaussian_poly_integral(la+lb,zet_a*zet_b/(zet_a+zet_b),r)/(pow(2*zet_a,1.5+la)*pow(2*zet_b,1.5+lb));
}
double prim_ovlp(double* ra,double* rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb)
{
   double result(0);

   double rab(sqrt(pow((ra[0]-rb[0]),2)+pow((ra[1]-rb[1]),2)+pow((ra[2]-rb[2]),2)));
   
   double thet( acos( (rb[2]-ra[2]) / rab ) );
   double phi( atan2( rb[1]-ra[1],rb[0]-ra[0] ) );

   for(unsigned int l=abs(la-lb);l!=la+lb+1;l++)
   {
      if((l+la+lb)%2!=0)
         continue;
      else
      {
         for(int m=-l;m<int(l)+1;m++)
         {
            if( (ma+mb+m) !=0 )
               continue;
            else
               result+=4*acos(-1)*pow(-1,(lb-la-l)/2)*prim_radial_ovlp(la,lb,zeta_a,zeta_b,rab)
               *rYlm(l,m,thet,phi)*three_Ylm_integ(la,lb,l,ma,mb,m);
         }
      }
   }
   return result;
}
double spherical_overlap_integral(double xa,double xb,double m,double p)
{
   return sqrt(acos(-1)/p)*exp(-m*fabs(xb-xa));
}
double obara_saika_ovlp(double xa,double xb,double zeta_a,double zeta_b,int la,int lb)
{

   if(la<0 || lb<0)
      return 0;

   double p(zeta_a+zeta_b); 
   double m(zeta_a*zeta_b/(zeta_a+zeta_b)); 
   double xp((zeta_a*xa+zeta_b*xb)/p);

   if(la==0 && lb==0)
      return spherical_overlap_integral(xa,xb,m,p);

   double xpa(xp-xa);
   double xpb(xp-xb);
   double Sam1b(obara_saika_ovlp(xa,xb,zeta_a,zeta_b,la-1,lb));
   double Sabm1(obara_saika_ovlp(xa,xb,zeta_a,zeta_b,la,lb-1));
   double Sam2b(obara_saika_ovlp(xa,xb,zeta_a,zeta_b,la-2,lb));
   double Sabm2(obara_saika_ovlp(xa,xb,zeta_a,zeta_b,la,lb-2));

   if(la >= lb)
      return  xpa*Sam1b + (1/(2*p))*( (la-1) * Sam2b + (lb-1) * Sabm1 );
   else
      return  xpb*Sabm1 + (1/(2*p))*( (la-1) * Sam1b + (lb-1) * Sabm2 );



}
