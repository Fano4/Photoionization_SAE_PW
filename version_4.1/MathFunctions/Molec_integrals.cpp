#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <arb.h>
#include <acb_hypgeom.h>
#include "mathfunctions.h"
//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_errno.h>

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


/*double bessel_gaussian_poly_integral(unsigned int l1,unsigned int l2,double m,double r)
//{
   //Computes the indefinite integral of k**(2+l1)*exp(-p*k*k)*j_l2(k*r)

//   std::cout<<"####"<<l1<<","<<l2<<","<<m<<","<<r<<std::endl;
   
   arb_t n;
   arb_init(n);
   arb_t a;
   arb_init(a);
   arb_t tmp;
   arb_init(tmp);

   arb_t thresh;
   arb_init(thresh);
   arb_set_d(thresh,1e-15);
   slong prec(256);
   
   
   arb_set_d(a,double(m/(r*r)));
   arb_set_d(n,(l1+l2)/2.);
   //using power series
   arb_t t;
   arb_init(t); 
   arb_t tt;
   arb_init(tt);
   arb_t tmp2;
   arb_init(tmp2);
   arb_t temp_a;
   arb_init(temp_a);
   arb_t temp_b;
   arb_init(temp_b);
   arb_t sum;
   arb_init(sum); 
   arb_t S;
   arb_init(S); 
   arb_t prefac;
   arb_init(prefac); 

   arb_set_d(S,1.); // S=1
   arb_set_d(tmp,1.); // S=1
   arb_set_d(prefac,pow(r,-l1-3)*sqrt(acos(-1))/(pow(2,l2+2)*pow(arf_get_d(arb_midref(a),ARF_RND_NEAR),1.5+arf_get_d(arb_midref(n),ARF_RND_NEAR)))); // prefactor=pow(r,-l1-3)*sqrt(acos(-1))/(pow(2,l2+2)*pow(a,1.5+n))

   arb_zero(sum); //sum=0
   arb_zero(t); //t=0

   std::cout<<std::endl<<"vvvv"<<std::endl;
   arb_print(prefac);
   std::cout<<std::endl;
   arb_print(a);
   std::cout<<std::endl;
   arb_print(n);
   std::cout<<std::endl;


   while( arf_cmpabs(arb_midref(tmp),arb_midref(thresh)) >= 0 )
   {
      arb_set_d(temp_a,1.5);
      arb_add(temp_a,temp_a,n,prec);
      arb_add(temp_a,temp_a,t,prec);

      arb_set_d(temp_b,1.5+l2);
      arb_add(temp_b,temp_b,t,prec);
      arb_set_d(S,1.-2*(int(arf_get_d(arb_midref(t),ARF_RND_NEAR))%2)); //pow(-1,t)
      arb_gamma(tmp,temp_a,prec);
      arb_mul(S,S,tmp,prec);
      arb_gamma(tmp,temp_b,prec);
      arb_div(S,S,tmp,prec);

      arb_set_d(tmp2,2);
      arb_mul(tmp,tmp2,t,prec);
      arb_pow(tmp,tmp2,tmp,prec);
      arb_div(S,S,tmp,prec);
      arb_pow(tmp,a,t,prec);
      arb_div(S,S,tmp,prec);
      arb_mul(S,S,prefac,prec);
      arb_set_d(tmp,1.);
      arb_add(tmp,t,tmp,prec);
      arb_gamma(tmp,tmp,prec);
      arb_div(S,S,tmp,prec);
      arb_add(sum,sum,S,prec);
      arb_set_d(tmp2,1);
      arb_add(t,t,tmp2,prec);

      arb_div(tmp,S,sum,prec);
      
      std::cout<<arf_get_d(arb_midref(t),ARF_RND_NEAR)<<"=>"<<"+++"<<arf_get_d(arb_midref(S),ARF_RND_NEAR)<<std::endl;
      if(arf_get_d(arb_midref(t),ARF_RND_NEAR) > 1000)
      {
         std::cout<<"No convergence for overlap integral. exit"<<std::endl;
         break;
      }

   }

   return arf_get_d(arb_midref(sum),ARF_RND_NEAR);
   
}*/
double prim_radial_ovlp(unsigned int la,unsigned int lb,unsigned int l,double zet_a,double zet_b,double r)
{

   //Special case if la+lb==l
/*   if(la+lb==l)
   {
      double m(zet_a*zet_b/(zet_a+zet_b));
      double a(1./(4.*m));
      return sqrt(acos(-1))*exp(-r*r/(4.*a))*pow(r,l)/(pow(2.,l+2)*pow(a,1.5+l)*pow(zet_a,1.5+la)*pow(zet_b,1.5+lb)); 
   }*/

   slong prec(128);

   arb_t ala,alb,al,azeta,azetb,ar; // Variables from input
   arb_init(ala);
   arb_init(alb);
   arb_init(al);
   arb_init(azeta);
   arb_init(azetb);
   arb_init(ar);
   arb_set_d(ala,int(la));
   arb_set_d(alb,int(lb));
   arb_set_d(al,int(l));
   arb_set_d(azeta,zet_a);
   arb_set_d(azetb,zet_b);
   arb_set_d(ar,r);

   arb_t prefac,tmp,tmp2,tmp3,tmp4,ares; //computation variables
   arb_init(prefac);
   arb_init(tmp);
   arb_init(tmp2);
   arb_init(tmp3);
   arb_init(tmp4);
   arb_init(ares);


   arb_add(tmp,alb,al,prec); // tmp= lb+l
   arb_sub(tmp,tmp,ala,prec); // tmp= lb+l-la
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (lb+l-la)/2

   arb_pow(prefac,azeta,tmp,prec); // prefac= pow(zet_a,(lb+l-la)/2)

   arb_add(tmp,ala,al,prec); // tmp= la+l
   arb_sub(tmp,tmp,alb,prec); // tmp= la+l-lb
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (la+l-lb)/2
   arb_pow(tmp,azetb,tmp,prec); // tmp= pow(zet_b,(la+l-lb)/2)

   arb_mul(prefac,prefac,tmp,prec); //prefac= pow(zet_a,(lb+l-la)/2) * pow(zet_b,(lb+l-la)/2)

   arb_add(tmp,ala,al,prec); // tmp= la+l
   arb_add(tmp,tmp,alb,prec); // tmp= la+l+lb
   arb_set_d(tmp2,3.);
   arb_add(tmp,tmp,tmp2,prec); // tmp= la+l+lb+3
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (la+l+lb+3)/2
   arb_add(tmp2,azeta,azetb,prec);
   arb_pow(tmp,tmp2,tmp,prec); // tmp= pow(zet_a+zet_b,(la+l+lb+3)/2)

   arb_div(prefac,prefac,tmp,prec); // prefac= pow(zet_a,(lb+l-la)/2) * pow(zet_b,(lb+l-la)/2) / pow(zet_a+zet_b,(la+l+lb+3)/2)

   arb_set_d(tmp,sqrt(acos(-1))/4);

   arb_mul(prefac,prefac,tmp,prec); // prefac= sqrt(pi)/4 * pow(zet_a,(lb+l-la)/2) * pow(zet_b,(lb+l-la)/2) / pow(zet_a+zet_b,(la+l+lb+3)/2)

   arb_pow(tmp,ar,al,prec); //r**l

   arb_mul(prefac,prefac,tmp,prec);


   arb_add(tmp,ala,al,prec); // tmp= la+l
   arb_add(tmp,tmp,alb,prec); // tmp= la+l+lb
   arb_set_d(tmp2,3.);
   arb_add(tmp,tmp,tmp2,prec); // tmp= la+l+lb+3
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (la+l+lb+3)/2
   arb_gamma(tmp,tmp,prec);

   arb_mul(prefac,prefac,tmp,prec); //prefac*=gamma((la+l+lb+3)/2)

   //Now setting up the hypergeometric

   arb_add(tmp,ala,al,prec); // tmp= la+l
   arb_add(tmp,tmp,alb,prec); // tmp= la+l+lb
   arb_set_d(tmp2,3.);
   arb_add(tmp,tmp,tmp2,prec); // tmp= la+l+lb+3
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (la+l+lb+3)/2

   arb_set_d(tmp2,1.5);
   arb_add(tmp2,tmp2,al,prec); // tmp2=1.5+l

   arb_add(tmp4,azeta,azetb,prec); // tmp4=zet_a+zet_b
   arb_mul(tmp3,azeta,azetb,prec); //tmp3=zet_a*zet_b
   arb_div(tmp3,tmp3,tmp4,prec); //tmp3=zet_a*zet_b / (zet_a+zet_b)
   arb_set_d(tmp4,-1.); 
   arb_mul(tmp3,tmp3,tmp4,prec); // tmp3=-zet_a*zet_b / (zet_a+zet_b)
   arb_mul(tmp4,ar,ar,prec); //tmp4=r**2
   arb_mul(tmp3,tmp3,tmp4,prec); // tmp3= - r**2 *zet_a*zet_b / (zet_a+zet_b)

   acb_t ctmp,ctmp2,ctmp3,cres;
   acb_init(ctmp);
   acb_init(ctmp2);
   acb_init(ctmp3);
   acb_init(cres);

   acb_set_arb(ctmp,tmp);
   acb_set_arb(ctmp2,tmp2);
   acb_set_arb(ctmp3,tmp3);
   acb_hypgeom_m(cres,ctmp,ctmp2,ctmp3,1,prec);

   acb_get_real(ares,cres);

   arb_mul(ares,ares,prefac,prec);

   return arf_get_d(arb_midref(ares),ARF_RND_NEAR);
}
double prim_ovlp(std::vector<double> ra,std::vector<double> rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb)
{

   double result(0);
   double thet;
   double phi;
   double temp(0);

   double rab(sqrt(pow((ra.at(0)-rb.at(0)),2.)+pow((ra.at(1)-rb.at(1)),2.)+pow((ra.at(2)-rb.at(2)),2.)));
   
   if(rab!=0)
   {
      thet=( acos( (rb.at(2)-ra.at(2)) / rab ) );
      if( (rb.at(0)-ra.at(0)) != 0)
         phi=(atan2( rb.at(1)-ra.at(1),rb.at(0)-ra.at(0) ) );
      else
         phi=0;
   }
   else
   {
      thet=(0);
      phi=(0);
   }

   //std::cout<<std::defaultfloat;
//   std::cout<<std::endl<<std::scientific<<"###"<<rab<<","<<thet<<","<<phi<<","<<zeta_a<<","<<zeta_b<<","<<la<<","<<lb<<","<<ma<<","<<mb<<std::endl<<"***";
   for(unsigned int l=abs(int(la-lb));l<=la+lb;l++)
   {
      if((l+la+lb)%2!=0)
         continue;
      else
      {
         // compute the factor independent of m
         temp=4*acos(-1)*pow(-1,(lb-la+l)/2)*prim_radial_ovlp(la,lb,l,zeta_a,zeta_b,rab);

         //Add the m=0 term
         result+=temp*rYlm(l,0,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,0)
               *three_azim_integ(ma,mb,0)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),0);

         //ass the other terms of the sum
         for(int m=1;m<=int(l);m++)
         {
            if( (ma+mb+m)%2!=0 )
               continue;
            else
               result+=2*temp*(rYlm(l,m,thet,phi)
               *prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,m)
               *three_azim_integ(ma,mb,m)
               *three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),m)
               +rYlm(l,-m,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,-m)
               *three_azim_integ(ma,mb,-m)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),-m));//*three_Ylm_integ(la,lb,l,ma,mb,m);
         }
      }
   }
//   std::cout<<" ==> "<<result<<std::endl;
   return result;
}

double ao_ovlp(std::vector<double> ra,std::vector<double> rb,std::vector<double> zet_a,std::vector<double> zet_b,std::vector<double> cont_coeff_a, std::vector<double> cont_coeff_b,unsigned int la,unsigned int lb,int ma,int mb)
{
   double result(0);
   for(unsigned int a=0;a!=cont_coeff_a.size();a++)
   {
      for(unsigned int b=0;b!=cont_coeff_b.size();b++)
      {
         result+=cont_coeff_a.at(a)*cont_coeff_b.at(b)*prim_ovlp(ra,rb,zet_a.at(a),zet_b.at(b),la,lb,ma,mb);
      }
   }
   return result;
}
/*
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
*/
