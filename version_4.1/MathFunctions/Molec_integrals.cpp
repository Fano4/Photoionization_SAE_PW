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
#include "utilities.h"

////////////////////////////////////////////////////////
//
//This computes the integrals denoted byt \mathcal{R} in the Latex document. It can be used for both computing overlap and dipole integrals
//It basically computes the radial part of the overlap integral between primitives
//
///////////////////////////////////////////////////////
double prim_radial_ovlp(unsigned int la,unsigned int lb,unsigned int l,double zet_a,double zet_b,double r)
{
   double m(zet_a*zet_b/(zet_a+zet_b));


   return gen_I_integ(la+lb,l,1./(4.*m),r)/(pow(2*zet_a,1.5+la)*pow(2*zet_b,1.5+lb));
   //return (la>=0 && lb>=0) ? gen_I_integ(la+lb,l,1./(4.*m),r)/(pow(2*zet_a,1.5+la)*pow(2*zet_b,1.5+lb)):0 ;
}
///////////////////////////////////////////////////////
//
//This computes the overlap between primitive GTOs. They include both the radial part and the angular parts and are the backend of
//the AO overlap integral computation. 
//
//////////////////////////////////////////////////////
double prim_ovlp(std::vector<double> ra,std::vector<double> rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb)
{

   double result(0);
   double thet;
   double phi;
   double temp(0);

   //The primitives are not normalized from their contraction coefficients. The normalization of each primitive is imposed here.

   double norma(sqrt(.5*tgamma(1.5+la)/pow(2*zeta_a,1.5+la)));
   double normb(sqrt(.5*tgamma(1.5+lb)/pow(2*zeta_b,1.5+lb)));

   //First compute the distance between both center. Zero distance meand that the overlap is betweeen basis functions located on the samne atom,
   //which leads to the simplest selection rule for the angular part

   double rab(sqrt(pow((ra.at(0)-rb.at(0)),2.)+pow((ra.at(1)-rb.at(1)),2.)+pow((ra.at(2)-rb.at(2)),2.)));
   
   if(rab!=0)
   {
      thet=( acos( (rb.at(2)-ra.at(2)) / rab ) );
      if( (rb.at(0)-ra.at(0)) != 0)
         phi=(atan2( rb.at(1)-ra.at(1),rb.at(0)-ra.at(0) ) );
      else
         phi=0;
   }
   else if(zeta_a == zeta_b)
      return bool(la==lb)*bool(ma==mb);
   else
   {
      return 0.5*tgamma(1.5+la)/(pow(zeta_a+zeta_b,1.5+la))*bool(la==lb)*bool(ma==mb)/(norma*normb); 
   }

   //If the distance is not zero, then compute the integral in reciprocal space, by expanding the phase factor in the basis of 
   //spherical harmonics. This leads to integrals over three spherical harmonics, and thus to selection rules.
   //This is an analytic result that involves wigner 3j symbols and gamma functions, but  no other special function
   result=0;
   for(unsigned int l=abs(int(la-lb));l<=la+lb;l++)
   {
      if((l+la+lb)%2!=0)
         continue;
      else
      {
         // compute the factor independent of m
         temp=4*acos(-1)*std::real(pow(std::complex<double>(0,-1),(la-lb-l)))*prim_radial_ovlp(la,lb,l,zeta_a,zeta_b,rab)/(norma*normb);

         //Add the m=0 term
         if( (ma+mb) % 2 == 0)
         {
            result+=temp*rYlm(l,0,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,0)
               *three_azim_integ(ma,mb,0)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),0);
         }

         //add the other terms of the sum
         for(int m=1;m<=int(l);m++)
         {
            if( (ma+mb+m)%2!=0 )
               continue;
            else
               result+=temp*(
               rYlm(l,m,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,m)
               *three_azim_integ(ma,mb,m)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),m)
               +rYlm(l,-m,thet,phi)*prefactor_rYlm(la,ma)*prefactor_rYlm(lb,mb)*prefactor_rYlm(l,-m)
               *three_azim_integ(ma,mb,-m)*three_ALP_J_integral(la,lb,l,abs(ma),abs(mb),m));//*three_Ylm_integ(la,lb,l,ma,mb,m);
         }
      }
   }
   return result;
}
///////////////////////////////////////////////////////
//
//This computes the overlap between atomic orbitals. The orbitals are GTOs written as contractions of primitives GTOs
//
//
///////////////////////////////////////////////////////

void ao_ovlp(std::vector<double> ra,std::vector<double> rb,std::vector<int> nuc_bas_fun_a,std::vector<int> nuc_bas_fun_b,std::vector<int> cont_num_a,std::vector<int> cont_num_b,std::vector<double> zet_a,std::vector<double> zet_b,std::vector<double> cont_coeff_a, std::vector<double> cont_coeff_b,std::vector<unsigned int> la,std::vector<unsigned int> lb,std::vector<int> ma,std::vector<int> mb,std::vector<double>* S)
{
   S->clear();
   double result(0);
   int counta(0);
   int countb(0);
   int mema(0);
   int memb(0);
   double temp(0);
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
/////////////////////////
//
//Computes the overlap between molecular orbitals by transforming the AO overlap matrix using the LCAO coefficients.
//
////////////////////////
void MO_ovlp(std::vector<double> S,std::vector<double> lcao_a,std::vector<double> lcao_b,std::vector<double>* MO_S)
{
   MO_S->clear();
   double* Ca=lcao_a.data();
   double* Cb=lcao_b.data();
   int basis_size(int(sqrt(S.size())));
   double* O=new double[basis_size*basis_size];

   int n_occ(int(lcao_a.size())/basis_size);
   double* res=new double [n_occ*n_occ];
   double* temp=new double [n_occ*basis_size];
   double* temp2=new double [n_occ*basis_size];

   //The AO necessitate a renormalization because the contraction coefficients in front of some primitives do not guarantee normalization.
   //This renormalization ensures the AO overlap matrix is the same as the one obtained in molpro.
   //
    for(int ao1=0;ao1!=basis_size;ao1++)
    {
       for(int ao2=0;ao2!=basis_size;ao2++)
       {
          O[ao1*basis_size+ao2]=S.at(ao1*basis_size+ao2)/sqrt(S.at(ao1*basis_size+ao1)*S.at(ao2*basis_size+ao2));
       }
    }

    //Matrix transformation
    transpose(Ca, temp2, n_occ, basis_size);
    matrix_product(temp, O, temp2, basis_size, basis_size, n_occ); 
    matrix_product(res, Cb, temp, n_occ,basis_size,n_occ);
    transpose(res,res, n_occ, n_occ);


    //record the result
    for(int i=0;i!=n_occ*n_occ;i++)
       MO_S->push_back(res[i]);

    delete [] temp2;
    delete [] temp;
    delete [] res;
    delete [] O;
}
void ES_ovlp(std::vector<double> CSF_S,int n_csf_a,int n_csf_b,std::vector<double> ci_vector_a,std::vector<double> ci_vector_b,int n_states_a,int n_states_b,std::vector<double>* ES_S)
{
   std::cout<<"Entering ES overlap routine with "<<n_csf_a<<" * "<<n_csf_b<<" = "<<CSF_S.size()<<" CSFs"<<std::endl;
   std::cout<<"There are "<<n_states_a<<" * "<<n_states_b<<" = "<<n_states_a*n_states_b<<" elements to compute"<<std::endl;
   ES_S->clear();
   double* Ca=ci_vector_a.data();
   double* Cb=ci_vector_b.data();
   double* O=CSF_S.data();
   double res[n_states_a*n_states_b];
   double temp[n_states_b*n_csf_a];
   double temp2[n_states_a*n_csf_a];


   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n_csf_a,n_states_b,n_csf_b,1,O,n_csf_b,Cb,n_states_b,0,&temp[0],n_states_b);
   cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,n_states_a,n_states_b,n_csf_a,1,Ca,n_states_a,&temp[0],n_states_b,0,&res[0],n_states_b);

    for(int i=0;i!=n_states_a*n_states_b;i++)
       ES_S->push_back(res[i]);

    std::cout<<" Exiting ES overlap routine"<<std::endl;
    return;
}
void matrix_product(double *C,double *A,double *B,const int dim1,const int dim2,const int dim3)
{
   //C=A*B
   std::cout<<" In matrix product with dimensions "<<dim1<<","<<dim2<<","<<dim3<<std::endl;
    double ntemp;
    for (int i=0; i<dim1; i++)
    {
        for (int j=0; j<dim3; j++)
        {
            ntemp=0;
            
            for (int k=0; k<dim2; k++)
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
double gen_I_integ(unsigned int l1,unsigned int l2,double zeta,double k)
{
   //This function compute the radial integral int_0^{inf} dr r^{2+l1}*exp(-d r^2)*j_{l2}(r * k)
   //For all practical purposes in this code, this function can only handle cases where l1+l2 is even
   //
   //error handling
   if(k*k/(4.*zeta) > 708.4)
      std::cout<<"Warning ! Risk of underflow in evaluating exponential function. x = "<<-k*k/(4.*zeta)<<" for "<<zeta<<" and "<<k<<std::endl;

   //computation 
   if(l1==l2)
      return sqrt(acos(-1)/2.)*pow(k,l1)*exp(-k*k/(4.*zeta))/pow(2*zeta,1.5+l1); 

   else if( abs(int(l1-l2)) % 2 == 0 && l1 > l2)
      return ((2*l2+3)/k)*gen_I_integ(l1-1,l2+1,zeta,k)-gen_I_integ(l1,l2+2,zeta,k);

   else if( abs(int(l1-l2)) % 2 == 0 && l1 < l2)
      return ((l1+l2+1)/k)*gen_I_integ(l1-1,l2-1,zeta,k)-(2*zeta/k)*gen_I_integ(l1+1,l2-1,zeta,k);

   else //if( abs(l1-l2) % 2 == 0 && l1 < l2)
      err_bad_indices_gen_I_integ(l1,l2);

   return 0;

}
//////////////////////////////////////////////////
//
//This computes the trtansition dipole between two primitives. The origin of the frame is determined by the coordinates of atom A and B. 
//The function accounts for the shift of atom A with respect to the origin
//
//////////////////////////////////////////////////
void prim_trdip(std::vector<double> ra,std::vector<double> rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb,std::vector<double>* trdip)
{
   std::vector<double> trdip_centered;
   double Sab;

   prim_trdip_centered(ra,rb,zeta_a,zeta_b,la,lb,ma,mb,&trdip_centered);
   Sab=prim_ovlp(ra,rb,zeta_a,zeta_b,la,lb,ma,mb);

   trdip->clear();

   trdip->push_back( trdip_centered.at(0) + ra.at(0) * Sab );
   trdip->push_back( trdip_centered.at(1) + ra.at(1) * Sab );
   trdip->push_back( trdip_centered.at(2) + ra.at(2) * Sab );

   return;
}
//////////////////////////////////////////////////
//
//This computes the transition dipole between two primitives. The dipole is computed in the frame cenetered on atom A.
//This exploits the same machinery and relations as the overlap.
//Because the selection rules are affected by the dipole operator, there are a few more terms in the sum. This implementation tries
//to deal with it in a simplistic way, but we should improve this
//
//////////////////////////////////////////////////
void prim_trdip_centered(std::vector<double> ra,std::vector<double> rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb,std::vector<double>* trdip)
{

   std::cout<<"Entering prim_trdip_centered with parameters ra=("<<ra.at(0)<<","<<ra.at(1)<<","<<ra.at(2)
      <<") ; rb=("<<ra.at(0)<<","<<ra.at(1)<<","<<ra.at(2)<<") ; "
      <<", zeta_a = "<<zeta_a<<" , zeta_b = "<<zeta_b
      <<", la = "<<la<<" , lb = "<<lb
      <<", ma = "<<ma<<" , mb = "<<mb<<std::endl;
   double tempx,tempy,tempz;
   double Bm1m1,Bm1p1,B0m1,B0p1,Bp1m1,Bp1p1;
   double Rlp1, Rlm1;
   double prefactor;
   double thet(0);
   double phi(0);
   std::vector<double> B_temp;

   //The primitives are not normalized from their contraction coefficients. The normalization of each primitive is imposed here.

   double norma(sqrt(.5*tgamma(1.5+la)/pow(2*zeta_a,1.5+la)));
   double normb(sqrt(.5*tgamma(1.5+lb)/pow(2*zeta_b,1.5+lb)));

   //First compute the distance between both center. Zero distance means that the dipole is betweeen basis functions located on the samne atom,
   //which leads to the simplest selection rule for the angular part

   double rab(sqrt(pow((ra.at(0)-rb.at(0)),2.)+pow((ra.at(1)-rb.at(1)),2.)+pow((ra.at(2)-rb.at(2)),2.)));

   B_coeff(la,ma,-1,&B_temp);
   Bm1m1=B_temp.at(0);
   Bm1p1=B_temp.at(1);
   B_coeff(la,ma,0,&B_temp);
   B0m1=B_temp.at(0);
   B0p1=B_temp.at(1);
   B_coeff(la,ma,1,&B_temp);
   Bp1m1=B_temp.at(0);
   Bp1p1=B_temp.at(1);


   if(rab!=0)
   {
      thet=( acos( (rb.at(2)-ra.at(2)) / rab ) );
      if( (rb.at(0)-ra.at(0)) != 0)
         phi=(atan2( rb.at(1)-ra.at(1),rb.at(0)-ra.at(0) ) );
   }
//   else if(zeta_a == zeta_b)
//      return bool(la==lb)*bool(ma==mb);
//   else
//   {
      //Special cases when on a single center???
//      return 0.5*tgamma(1.5+la)/(pow(zeta_a+zeta_b,1.5+la))*bool(la==lb)*bool(ma==mb)/(norma*normb); 
//   }

   //If the distance is not zero, then compute the integral in reciprocal space, by expanding the phase factor in the basis of 
   //spherical harmonics. This leads to integrals over three spherical harmonics, and thus to selection rules.
   //This is an analytic result that involves wigner 3j symbols and gamma functions, but  no other special function

   tempx=0;
   tempy=0;
   tempz=0;

   for(unsigned int l=abs(int(la-lb+1));l<=la+lb+1;l++)
   {
      std::cout<<" l = "<<l<<std::endl;
      if((l+la+lb+1)%2==0)
      {
         Rlm1=prim_radial_ovlp(la-1,lb,l,zeta_a,zeta_b,rab);
         Rlp1=prim_radial_ovlp(la+1,lb,l,zeta_a,zeta_b,rab);
         // compute the factor independent of m
         if(ma==0)
         {
            prefactor = pow(4*acos(-1)/3.,1.5) * 3 * std::real(pow(std::complex<double>(0,-1),(la-lb+l+1)))/(norma*normb);

            tempx += prefactor * ( B0m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,1,mb) //0 if la+lb+l+1 is odd
                                 + B0p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,1,mb) //0 if la+lb+l+1 is odd
                                 - ( ( 2 * la + 1) ) * ( B0m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,1,mb) ) //0 if la+lb+l+1 is odd
                                 );

            tempy += prefactor * ( B0m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,-1,mb) 
                                 + B0p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,-1,mb)
                                 - ( ( 2 * la + 1) ) * ( B0m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,-1,mb) )
                                 );

            tempz += prefactor * ( B0m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,0,mb) 
                                 + B0p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,0,mb) 
                                 - ( ( 2 * la + 1) ) * ( B0m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,0,mb) )
                                 );
         }
         else
         {
            prefactor = pow(4*acos(-1)/3.,1.5) * (3. / sqrt(2.) ) * std::real(pow(std::complex<double>(0,-1),(la-lb+l+1)))/(norma*normb);

            tempx += prefactor * ( Bm1m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,(ma/abs(ma))*(abs(ma)-1),mb) 
                                 + Bp1m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,(ma/abs(ma))*(abs(ma)+1),mb)
                                 + Bm1p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,(ma/abs(ma))*(abs(ma)-1),mb)
                                 + Bp1p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,(ma/abs(ma))*(abs(ma)+1),mb)
                                 - ( ( 2 * la + 1 ) ) * ( Bp1m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,(ma/abs(ma))*(abs(ma)+1),mb)
                                                             + Bm1m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,(ma/abs(ma))*(abs(ma)-1),mb))
                                 );

            tempy += -(ma/abs(ma)) * prefactor * ( - Bm1m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,-(ma/abs(ma))*(abs(ma)-1),mb) 
                                                 + Bp1m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,-(ma/abs(ma))*(abs(ma)+1),mb)
                                                 - Bm1p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,-(ma/abs(ma))*(abs(ma)-1),mb)
                                                 + Bp1p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,-(ma/abs(ma))*(abs(ma)+1),mb)
                                                 - ( ( 2 * la + 1 ) ) * ( Bp1m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,-(ma/abs(ma))*(abs(ma)+1),mb)
                                                                             - Bm1m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,-(ma/abs(ma))*(abs(ma)-1),mb))
                                                 );

            tempz += sqrt(2) * prefactor * ( B0m1 * Rlp1 * CY_m_sum(thet,phi,la-1,lb,l,ma,mb) 
                                           + B0p1 * Rlp1 * CY_m_sum(thet,phi,la+1,lb,l,ma,mb) 
                                           - ( ( 2 * la + 1) ) * ( B0m1 * Rlm1 * CY_m_sum(thet,phi,la-1,lb,l,ma,mb) )
                                           );
         }
      }
   }

   trdip->clear();
   trdip->push_back(tempx);
   trdip->push_back(tempy);
   trdip->push_back(tempz);

   return;

}
