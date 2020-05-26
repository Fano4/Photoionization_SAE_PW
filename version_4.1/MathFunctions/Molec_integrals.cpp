#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <complex>
#include <arb.h>
#include <arb_hypgeom.h>
#include <mkl.h>
#include <omp.h>
#include "mathfunctions.h"
//#include <gsl/gsl_sf_hyperg.h>
//#include <gsl/gsl_errno.h>

/*
bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff)
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
}
*/

double prim_radial_ovlp(unsigned int la,unsigned int lb,unsigned int l,double zet_a,double zet_b,double r)
{

   //Special case if la+lb==l
   if(la+lb==l)
   {
      double m(zet_a*zet_b/(zet_a+zet_b));
      return sqrt(acos(-1))*exp(-m*r*r)*pow(r,l)*pow(m,1.5+l)*pow(2.,l+1)/(pow(2*zet_a,1.5+la)*pow(2*zet_b,1.5+lb)); 
   }

   slong prec(256);

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

   arb_mul(prefac,prefac,tmp,prec); //prefac= pow(zet_a,(la+l-lb)/2) * pow(zet_b,(la+l-lb)/2)

   arb_add(tmp,ala,al,prec); // tmp= la+l
   arb_add(tmp,tmp,alb,prec); // tmp= la+l+lb
   arb_set_d(tmp2,3.);
   arb_add(tmp,tmp,tmp2,prec); // tmp= la+l+lb+3
   arb_set_d(tmp2,2.);
   arb_div(tmp,tmp,tmp2,prec); // tmp= (la+l+lb+3)/2
   arb_add(tmp2,azeta,azetb,prec);
   arb_pow(tmp,tmp2,tmp,prec); // tmp= pow(zet_a+zet_b,(la+l+lb+3)/2)

   arb_div(prefac,prefac,tmp,prec); // prefac= pow(zet_a,(la+l-lb)/2) * pow(zet_b,(la+l-lb)/2) / pow(zet_a+zet_b,(la+l+lb+3)/2)

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

   arb_hypgeom_m(ares,tmp,tmp2,tmp3,1,prec);
   /*
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
   */

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
      return 0.5*tgamma(1.5+la)/(pow(zeta_a+zeta_b,1.5+la))*bool(la==lb)*bool(ma==mb); 
   }

   for(unsigned int l=abs(int(la-lb));l<=la+lb;l++)
   {
      if((l+la+lb)%2!=0)
         continue;
      else
      {
         // compute the factor independent of m
         temp=4*acos(-1)*std::real(pow(std::complex<double>(0,-1),(la-lb-l)))*prim_radial_ovlp(la,lb,l,zeta_a,zeta_b,rab);

         //Add the m=0 term
         result+=temp*rYlm(l,0,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,0)
               *three_azim_integ(ma,mb,0)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),0);

         //add the other terms of the sum
         for(int m=1;m<=int(l);m++)
         {
            if( (ma+mb+m)%2!=0 )
               continue;
            else
               result+=2*temp*(
               rYlm(l,m,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,m)
               *three_azim_integ(ma,mb,m)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),m)
               +rYlm(l,-m,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,-m)
               *three_azim_integ(ma,mb,-m)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),-m));//*three_Ylm_integ(la,lb,l,ma,mb,m);
         }
      }
   }
   return result;
}

void ao_ovlp(std::vector<double> ra,std::vector<double> rb,std::vector<int> nuc_bas_fun_a,std::vector<int> nuc_bas_fun_b,std::vector<int> cont_num_a,std::vector<int> cont_num_b,std::vector<double> zet_a,std::vector<double> zet_b,std::vector<double> cont_coeff_a, std::vector<double> cont_coeff_b,std::vector<unsigned int> la,std::vector<unsigned int> lb,std::vector<int> ma,std::vector<int> mb,std::vector<double>* S)
{
   double result(0);
   int counta(0);
   int countb(0);
   int mema(0);
   int memb(0);
   std::vector<double> rao1;
   std::vector<double> rao2;
   for(unsigned int ao1=0;ao1<cont_num_a.size();ao1++)
   {

      rao1.clear();
      rao1.push_back(ra.at(3*nuc_bas_fun_a.at(ao1)));
      rao1.push_back(ra.at(3*nuc_bas_fun_a.at(ao1)+1));
      rao1.push_back(ra.at(3*nuc_bas_fun_a.at(ao1)+2));

      countb=0;
      memb=0;
      for(unsigned int ao2=0;ao2<cont_num_b.size();ao2++)
      {

         rao2.clear();
         rao2.push_back(rb.at(3*nuc_bas_fun_b.at(ao2)));
         rao2.push_back(rb.at(3*nuc_bas_fun_b.at(ao2)+1));
         rao2.push_back(rb.at(3*nuc_bas_fun_b.at(ao2)+2));
         //std::cout<<"AO1 "<<ao1<<" => nucleus "<<nuc_bas_fun_a.at(ao1)<<std::endl;
         //std::cout<<rao1.at(0)<<" , "<<rao1.at(1)<<" , "<<rao1.at(2)<<std::endl<<"+++"<<std::endl;
         //std::cout<<"AO2 "<<ao2<<" => nucleus "<<nuc_bas_fun_b.at(ao2)<<std::endl;
         //std::cout<<rao2.at(0)<<" , "<<rao2.at(1)<<" , "<<rao2.at(2)<<std::endl<<"+++"<<std::endl;

         result=0;
         counta=mema;
         for(unsigned int a=0;a<cont_num_a.at(ao1);a++)
         {
            countb=memb;
            for(unsigned int b=0;b<cont_num_b.at(ao2);b++)
            {
               result+=cont_coeff_a.at(counta)*cont_coeff_b.at(countb)*prim_ovlp(rao1,rao2,zet_a.at(counta),zet_b.at(countb),la.at(ao1),lb.at(ao2),ma.at(ao1),mb.at(ao2));
               countb++;
            }
            counta++;
         }
         memb=countb;
         S->push_back(result);
      }
      mema=counta;
   }
}
void MO_ovlp(std::vector<double> S,std::vector<double> lcao_a,std::vector<double> lcao_b,std::vector<double>* MO_S)
{
   double* Ca=lcao_a.data();
   double* Cb=lcao_b.data();
   double* O=S.data();
   int basis_size(int(sqrt(S.size())));
   int n_occ(int(lcao_a.size())/basis_size);
   double* res=new double [n_occ*n_occ];
   double* temp=new double [n_occ*basis_size];
   double* temp2=new double [n_occ*basis_size];
//   double temp(0);

/*
   for(int i=0;i!=n_occ;i++)
   {
      for(int j=0;j!=n_occ;j++)
      {
         temp=0;
         for(int r=0;r!=basis_size;r++)
         {
            for(int s=0;s!=basis_size;s++)
            {
               temp+=S.at(r*basis_size+s)*lcao_a.at(i*basis_size+r)*lcao_b.at(j*basis_size+s);
            }
         }
         MO_S->push_back(temp);
      }
   }*/
//   std::cout<<"+++"<<lcao_a.size()<<" / "<<basis_size<<" = "<<n_occ<<std::endl;


    transpose(Ca, temp2, n_occ, basis_size);
    matrix_product(temp, O, temp2, basis_size, basis_size, n_occ); 
    matrix_product(res, Cb, temp, n_occ,basis_size,n_occ);
    transpose(res,res, n_occ, n_occ);

    for(int i=0;i!=n_occ*n_occ;i++)
       MO_S->push_back(res[i]);

    delete [] temp2;
    delete [] temp;
    delete [] res;
}
void ES_ovlp(std::vector<double> CSF_S,int n_csf_a,int n_csf_b,std::vector<double> ci_vector_a,std::vector<double> ci_vector_b,int n_states_a,int n_states_b,std::vector<double>* ES_S)
{
   double* Ca=ci_vector_a.data();
   double* Cb=ci_vector_b.data();
   double* O=CSF_S.data();
   double* res=new double [n_states_a*n_states_b];
   double* temp=new double [n_states_b*n_csf_b];
   double* temp2=new double [n_states_a*n_csf_a];

    transpose(Cb, temp2,n_csf_a,n_states_a);
    matrix_product(temp, O, Ca, n_csf_b, n_csf_a, n_states_a); 
    matrix_product(res, temp2, temp, n_states_b, n_csf_b, n_states_a);
    transpose(res,res, n_states_b, n_states_a);

    for(int i=0;i!=n_states_a*n_states_b;i++)
       ES_S->push_back(res[i]);

    delete [] temp2;
    delete [] temp;
    delete [] res;
}
void matrix_product(double *C,double *A,double *B,int dim1,int dim2,int dim3)
{
   //C=A*B
    double ntemp;
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim3; j++)
        {
            ntemp=0;
            
            for (int k=0; k!=dim2; k++)
            {
                ntemp+=A[i*dim2+k]*B[k*dim3+j];
            }
            C[i*dim3+j]=ntemp;
        }
    }
    
}
void transpose(double *A,double *B, int dim1, int dim2)
{
    //B=trans(A)
   if(A!=B)
   {
     for (int i=0; i!=dim1; i++)
     {
        for (int j=0; j!=dim2; j++)
        {
            B[j*dim1+i]=A[i*dim2+j];
        }
     }
   }
   else
   {
      double* temp=new double[dim1*dim2];
      for (int i=0; i!=dim1; i++)
      {
         for (int j=0; j!=dim2; j++)
         {
             temp[j*dim1+i]=A[i*dim2+j];
         }
      }
      for(int i=0;i!=dim1*dim2;i++)
         A[i]=temp[i];

      delete [] temp;
    }
}
double determinant(double *A,int dim)
{
    double det_val(1);
    short int sign(1);
    double *B=new double[dim*dim];
    int *ipiv=new int[dim];
    int n(0);
    for (int i=0; i!=dim; i++)
    {
        for (int j=0; j!=dim; j++)
        {
            B[i*dim+j]=A[i*dim+j];
        }
    }

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR,dim,dim,B,dim,ipiv);

    for(int i=0;i!=dim;i++)
    {
       if(i+1!=ipiv[i])
          sign*=-1;

       det_val*=B[i*dim+i];
    }
    delete [] ipiv;
    delete [] B;
    return sign*det_val;
}
