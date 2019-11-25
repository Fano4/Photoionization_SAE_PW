#include "Computation.hpp"
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
/*
std::complex<double> MO_Fourier_transform_grad( int mo_index,int comp, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size)
{
//   std::cout<<"MO_FT_GRAD :"<<mo_index<<","<<comp<<","<<k<<","<<thet<<","<<phi<<","<<basis_size<<std::endl;
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      if(MO_neut_basis_coeff[mo_index*basis_size+i] != 0)
      {
          value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT_grad(i,comp,k,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,angular_mom_numbers[i]);
//          std::cout<<" => "<<MO_neut_basis_coeff[mo_index*basis_size+i]<<","<<nucl_basis_func[i]-1<<": "<<value<<std::endl;
      }
   }
//      exit(EXIT_SUCCESS);
//   std::cout<<value<<"=="<<std::endl;
   if(real(value) == NAN || imag(value) == NAN)
   {
      std::cout<<"FATAL ERROR. UNDEFINED COUPLING VALUE"<<std::endl<<value<<std::endl;
      exit(EXIT_FAILURE);
   }
   return value;
}
std::complex<double> AO_FT_grad(int ao_index,int comp,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers,unsigned long long int* fact_memo)
{
   std::complex<double> value(0,0);
   switch(comp)
   {
      case 0:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
//   std::cout<<k<<","<<thet<<","<<phi<<","<<nucl_spher_pos[0]<<","<<contraction_zeta[ao_index][i]<<","<<angular_mom_numbers[0]<<"-"<<angular_mom_numbers[1]<<" ::: "<<contraction_coeff[ao_index][i]<<std::endl;
          if(contraction_coeff[ao_index][i]!=0)
             value+=contraction_coeff[ao_index][i]*contraction_FT_grad_k(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers,fact_memo);
//         std::cout<<"k => "<<contraction_FT_grad_k(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers)<<std::endl;
       }
       break;
      case 1:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
//   std::cout<<k<<","<<thet<<","<<phi<<","<<nucl_spher_pos[0]<<","<<contraction_zeta[ao_index][i]<<","<<angular_mom_numbers[0]<<"-"<<angular_mom_numbers[1]<<" ::: "<<contraction_coeff[ao_index][i]<<std::endl;
          if(contraction_coeff[ao_index][i]!=0)
            value+=contraction_coeff[ao_index][i]*contraction_FT_grad_thet(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers,fact_memo);
//          std::cout<<"t => "<<contraction_FT_grad_thet(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers)<<std::endl;
       }
       break;
      case 2:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
//   std::cout<<k<<","<<thet<<","<<phi<<","<<nucl_spher_pos[0]<<","<<contraction_zeta[ao_index][i]<<","<<angular_mom_numbers[0]<<"-"<<angular_mom_numbers[1]<<" ::: "<<contraction_coeff[ao_index][i]<<std::endl;
          if(contraction_coeff[ao_index][i]!=0)
             value+=contraction_coeff[ao_index][i]*contraction_FT_grad_phi(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers,fact_memo);
//         std::cout<<"f => "<<contraction_FT_grad_phi(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers)<<std::endl;
       }
       break;
      default:
       std::cout<<"invalid component of grad vector required"<<std::endl;
       exit(EXIT_FAILURE);
       break;
   }
//   std::cout<<value<<"=="<<std::endl;
   return value;

}
std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,int **angular_mom_numbers,double *MO_neut_basis_coeff,int basis_size,unsigned long long int* fact_memo)
{
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      if(MO_neut_basis_coeff[mo_index*basis_size+i] != 0)
      {
         value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT(i,k,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,angular_mom_numbers[i],fact_memo);
      }
//      std::cout<<nucl_basis_func[i]-1<<"!!"<<nucl_spher_pos[nucl_basis_func[i]-1][0]<<","<<nucl_spher_pos[nucl_basis_func[i]-1][1]<<","<<nucl_spher_pos[nucl_basis_func[i]-1][2]<<std::endl;
   }
//   exit(EXIT_SUCCESS);
   return value;
}
std::complex<double> AO_FT(int ao_index,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,int* angular_mom_numbers,unsigned long long int* fact_memo)
{
   std::complex<double> value(0);
   for(int i=0;i!=contraction_number[ao_index];i++)
   {
      if(contraction_coeff[ao_index][i] != 0)
      {
         value+=contraction_coeff[ao_index][i]*contraction_FT(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],angular_mom_numbers,fact_memo);
      }
   }
   return value;

}
*/
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,double* lnfact_memo)
{
   using namespace std;

   int l(angular_mom_numbers[0]);
   int ml(angular_mom_numbers[1]);
   std::complex<double> value;
   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(std::complex<double>(0,-1)*k*nucl_spher_pos[0]*f));

   if(ml<0)
   {
   value=
      (sqrt(2)*associated_legendre(l,-ml,cos(thet))*sin(-ml*phi)) // Real spherical harmonics
      *((pow(std::complex<double>(0,-k),l)*exp(-k*k/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
   }
   else if(ml>0)
   {
   value=
      (sqrt(2)*associated_legendre(l,ml,cos(thet))*cos(ml*phi)) // Real spherical harmonics
      *((pow(std::complex<double>(0,-k),l)*exp(-k*k/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
      //std::cout<<thet<<";"<<phi<<" : probe! => "<<value<<std::endl;
   }
   else
   {
   value=
      (associated_legendre(l,ml,cos(thet))) // Real spherical harmonics
      *((pow(std::complex<double>(0,-k),l)*exp(-k*k/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
   }

   return std::conj(value);

}
/*
std::complex<double> contraction_FT_grad_k( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,unsigned long long int* fact_memo)
{
   using namespace std;

   int l(angular_mom_numbers[0]);
   int ml(angular_mom_numbers[1]);
   std::complex<double> value;
   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));
//   double *legendre_val=new double[gsl_sf_legendre_array_n(l)];
//   gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,l,cos(thet),-1,legendre_val);

//   std::cout<<"contraction_grad_k "<<nucl_spher_pos[0]<<";"<<contraction_zeta<<";"<<angular_mom_numbers[0]<<"-"<<angular_mom_numbers[1]<<std::endl;

   // ⎷(2) * P( l , |m| ) * (cos( m * phi )) * DK(k)/dk * exp( - I k r0 )  m > 0
   // ⎷(2) * P( l , |m| ) * (sin( |m| * phi )) * DK(k)/dk * exp( - I k r0 )  m < 0 


   if(ml!=0)
   {
//      std::cout<<"ml = "<<ml<<std::endl;
//
      value=
            sqrt(2)*
            (
 //             legendre_val[gsl_sf_legendre_array_index(l,fabs(ml))]
                associated_legendre(l,fabs(ml),cos(thet),fact_memo)
            ) // Associated Legendre
 
           *(
              bool(ml>0) * cos(ml * phi) + bool(ml<0) * sin(fabs(ml) * phi)
            ) // ml sign dependent phi part

            *(
               double(l) - k * ( k / (2*contraction_zeta) + std::complex<double>(0,1) * nucl_spher_pos[0] * f ) 
             )
            
            *(
               pow(std::complex<double>(0,-1),l) * pow(k,l-1) * exp( -k*k / (4*contraction_zeta) ) / pow(2*contraction_zeta,1.5+l)
             )
             
            *phase_factor // phase factor
      ;
//      std::cout<<"ddk : "<<value<<std::endl;
   }
   else
   {
//      std::cout<<"ml = "<<ml<<std::endl;
      value=
            (
              //legendre_val[gsl_sf_legendre_array_index(l,fabs(ml))]
                associated_legendre(l,fabs(ml),cos(thet),fact_memo)
            ) // Real spherical harmonics
      
            *(
               std::complex<double>(l,0) - k * ( k / (2*contraction_zeta) + std::complex<double>(0,1) * nucl_spher_pos[0] * f)
             )

            *(
               pow(std::complex<double>(0,-1),l) * pow(k,l-1) * exp( -k*k / (4*contraction_zeta) ) / pow(2*contraction_zeta,1.5+l)
             )
      
            *phase_factor // phase factor
      ;
//      std::cout<<"ddk : "<<value<<std::endl;
   }
   / *
   if(ml<0)
   {
   value=
      (pow(-1,ml)*sqrt(2)*(factorial(l-ml)/factorial(l+ml))*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*sin(-ml*phi)) // Real spherical harmonics
      *(std::complex<double>(l,0)-k*(k/(2*contraction_zeta)+std::complex<double>(0,1)*nucl_spher_pos[0]*f))*(pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta))/pow(2*contraction_zeta,1.5+l))
      *phase_factor // phase factor
      ;
   }
   else if(angular_mom_numbers[1]>0)
   {
   value=
      (pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*cos(ml*phi)) // Real spherical harmonics
      *(std::complex<double>(l,0)-k*(k/(2*contraction_zeta)+std::complex<double>(0,1)*nucl_spher_pos[0]*f))*(pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta))/pow(2*contraction_zeta,1.5+l))
      *phase_factor // phase factor
      ;
   }
   else
   {
   value=
      (gsl_sf_legendre_sphPlm(l,0,cos(thet))) // Real spherical harmonics
      *(std::complex<double>(l,0)-k*(k/(2*contraction_zeta)+std::complex<double>(0,1)*nucl_spher_pos[0]*f))*(pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta))/pow(2*contraction_zeta,1.5+l))
      *phase_factor // phase factor
      ;
   }
   * /
//   std::cout<<"dk:"<<l<<","<<ml<<","<<thet<<","<<phi<<","<<contraction_zeta<<","<<value<<std::endl;

//   delete [] legendre_val;
   return std::conj(value);

}
std::complex<double> contraction_FT_grad_thet( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,unsigned long long int* fact_memo)
{
   using namespace std;

   int l(angular_mom_numbers[0]);
   int ml(angular_mom_numbers[1]);
   std::complex<double> value;
   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdt(cos(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])-sin(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));
 //  double *legendre_val=new double[gsl_sf_legendre_array_n(l)];
 //  double *legendre_der_val=new double[gsl_sf_legendre_array_n(l)];

//   exit(EXIT_SUCCESS);

//   if(thet!=0 && thet!=acos(-1))
//   {
//      gsl_sf_legendre_deriv_alt_array_e(GSL_SF_LEGENDRE_SPHARM,l,cos(thet),-1,legendre_val,legendre_der_val);
//      std::cout<<l<<","<<ml<<" => "<<sqrt(3/(4*acos(-1)))*sin(thet)<<" | "<<legendre_val[gsl_sf_legendre_array_index(l,fabs(ml))]<<" d/dt = "<<legendre_der_val[gsl_sf_legendre_array_index(l,fabs(ml))]<<" - thet = "<<thet<<std::endl;
//   }
//   else
//   {
//      gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,l,cos(thet),-1,legendre_val);
//      for(int i=0;i!=gsl_sf_legendre_array_n(l);i++)
//         legendre_der_val[i]=0;
//   }
   
   // ⎷(2) * P( l , |m| ) * (cos( m * phi )) * DK(k)/dk * exp( - I * k * r0 * f)  m > 0
   // ⎷(2) * P( l , |m| ) * (sin( |m| * phi )) * DK(k)/dk * exp( - I * k * r0 * f)  m < 0


   if(ml!=0)
   {

//      std::cout<<"ml = "<<ml<<std::endl;
      value=
         sqrt(2)*
         (

///            std::complex<double>(0,-1) * k * nucl_spher_pos[0] * dfdt * legendre_val[gsl_sf_legendre_array_index(l,fabs(ml))] //Associated Legendre polynomial
            std::complex<double>(0,-1) * k * nucl_spher_pos[0] * dfdt * associated_legendre(l,fabs(ml),cos(thet),fact_memo) //Associated Legendre polynomial

             + associated_legendre_der(l,fabs(ml),cos(thet),fact_memo)
//           +  legendre_der_val[gsl_sf_legendre_array_index(l,fabs(ml))] 
         ) // Derivative of the associated Legendre polynomial. Theta dependent parts
 
         *( 
               bool(ml>0) * cos(ml * phi) + bool(ml<0) * sin(fabs(ml) * phi) 
          ) // ml sign dependent phi part

         *( 
               (pow(std::complex<double>(0,-1),l) * pow(k,l-1) * exp( -k*k / (4*contraction_zeta))) / (pow(2*contraction_zeta,1.5+l)) 
          ) // Radial part

          *phase_factor // phase factor
         ;
//      std::cout<<"ddt : "<<value<<std::endl;
   }
   else
   {
      //std::cout<<"ml = "<<ml<<std::endl;
      value=
         (std::complex<double>(0,-1) * k * nucl_spher_pos[0] * dfdt * associated_legendre(l,ml,cos(thet,fact_memo)) //Associated Legendre polynomial
            +associated_legendre_der(l,ml,cos(thet),fact_memo)) // Derivative of the associated Legendre polynomial. Theta dependent parts
         *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-k*k/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
         *phase_factor // phase factor
         ;
      //std::cout<<"ddt : "<<value<<std::endl;
   }
//std::cout<<"dt:"<<l<<","<<ml<<","<<thet<<","<<phi<<","<<contraction_zeta<<","<<value<<std::endl;
   / *
   if(ml<0)
   {
//      std::cout<<"probe 1 t : l = "<<l<<" ; ml = "<<ml<<std::endl;
      if(thet!=0)
      {
         value=
         ( std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(pow(-1,ml)*sqrt(2)*(factorial(l-ml)/factorial(l+ml))*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*sin(-ml*phi)))
//         +(pow(-1,ml)*sqrt(2)*(factorial(l-ml)/factorial(l+ml))*(sin(-ml*phi)/sin(thet))*(cos(thet)*l*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))-(l-ml)*sqrt((2*l+1)*(l+ml)/((2*l-1)*(l-ml)))*gsl_sf_legendre_sphPlm(l-1,bool(l+ml>0)*(-ml),cos(thet))*bool(l+ml>0))))// Real spherical harmonics
         *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
         *phase_factor // phase factor
         ;
      }
      else
      {
         value=
         ( std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(pow(-1,ml)*sqrt(2)*(factorial(l-ml)/factorial(l+ml))*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*sin(-ml*phi)))
         *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
         *phase_factor // phase factor
         ;
      }
   }
   else if(ml>0)
   {
//      std::cout<<"probe 2 t : l = "<<l<<" ; ml = "<<ml<<std::endl;
      if(thet!=0)
      {
         value=
           ( std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*cos(ml*phi))
           +(pow(-1,ml)*sqrt(2)*(cos(ml*phi)/sin(thet))*(cos(thet)*l*gsl_sf_legendre_sphPlm(l,ml,cos(thet))-(l+ml)*sqrt((2*l+1)*(l-ml)/((2*l-1)*(l+ml)))*gsl_sf_legendre_sphPlm(l-1,bool(l-ml>0)*(ml),cos(thet))*(bool(l-ml>0)))))// Real spherical harmonics
           *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
           *phase_factor // phase factor
           ;
      }
      else
      {
         value=
           ( std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*cos(ml*phi)))
           *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
           *phase_factor // phase factor
           ;
      }
   }
   else
   {
//      std::cout<<"probe 3 t : l = "<<l<<" ; ml = "<<ml<<std::endl;
      if(thet!=0)
      {
         value=
          (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*gsl_sf_legendre_sphPlm(l,0,cos(thet))*dfdt
          + (l/sin(thet))*(cos(thet)*gsl_sf_legendre_sphPlm(l,0,cos(thet))+gsl_sf_legendre_sphPlm(bool(l-1>=0)*(l-1),0,cos(thet))))// Real spherical harmonics
          *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
          *phase_factor // phase factor
          ;
      }
      else
      {
         value=
          (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*sqrt(2)*gsl_sf_legendre_sphPlm(l,0,cos(thet))*dfdt)
          *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
          *phase_factor // phase factor
          ;
      }
   }
* /
//   delete [] legendre_val;
//   delete [] legendre_der_val;
   return std::conj(value);
}
std::complex<double> contraction_FT_grad_phi( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,int* angular_mom_numbers,unsigned long long int* fact_memo)
{
   using namespace std;

   int l(angular_mom_numbers[0]);
   int ml(angular_mom_numbers[1]);
   std::complex<double> value;
   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdf(-sin(thet)*sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
   complex<double> dfdf2(-sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));

//      std::cout<<k<<" , "<<thet<<" , "<<phi<<" ; "<<contraction_zeta<<" , "<<l<<" , "<<ml<<","<<nucl_spher_pos[0]<<","<<nucl_spher_pos[1]<<","<<nucl_spher_pos[2]<<std::endl;
   if(ml!=0)
   {
      value=
         sqrt(2)*
         ( 
            (1/(2*fabs(ml)))*sqrt(((2*l+1)*factorial(l-fabs(ml)))/((4*acos(-1))*(factorial(l+fabs(ml))))) 
            *(
                 associated_legendre_nonorm(l+1,fabs(ml)+1,cos(thet),fact_memo) + (l-fabs(ml)+1) * (l-fabs(ml)+2) * associated_legendre_nonorm(l+1,fabs(ml)-1,cos(thet),fact_memo)
  //             legendre_val[gsl_sf_legendre_array_index(l+1,fabs(ml)+1)] + (l-fabs(ml)+1) * (l-fabs(ml)+2) * legendre_val[gsl_sf_legendre_array_index(l+1,fabs(ml)-1)]
             )

            *(

               fabs(ml) * bool(ml>0) * sin(ml * phi) - fabs(ml) * bool(ml<0) * cos(ml * phi)

               -std::complex<double>(0,-1) * k * nucl_spher_pos[0] * dfdf // Real spherical harmonics 
              
               *(
                  bool(ml>0) * cos(ml * phi) + bool(ml<0) * sin(fabs(ml) * phi)
                )
             )
          ) // ml sign dependent phi part

         *(
             pow(std::complex<double>(0,-1),l) * pow(k,l-1) * exp(-k*k / (4*contraction_zeta))  / (pow(2*contraction_zeta,1.5+l))
          ) // Radial part

         *phase_factor // phase factor
         ;
      
   }
   else
   {
//      std::cout<<"ml = "<<ml<<std::endl;
      value=
            (
               (
                   std::complex<double>(0,-1) * k *nucl_spher_pos[0] * associated_legendre(l,ml,cos(thet),fact_memo)
               ) // Real spherical harmonics
              *(
                 (pow(std::complex<double>(0,-1),l) * pow(k,l-1) * exp(-k*k / (4*contraction_zeta)) ) / (pow( 2 * contraction_zeta,1.5 + l)) 
               )
            ) // Radial part
            *phase_factor // phase factor
         ;
//      std::cout<<"ddf : "<<value<<std::endl;
   }
//   std::cout<<"df:"<<l<<","<<ml<<","<<thet<<","<<phi<<","<<contraction_zeta<<","<<value<<std::endl;
/ *
   if(ml<0)
   {
//      std::cout<<"probe 1 f : l = "<<l<<" ; ml = "<<ml<<std::endl;
//      std::cout<<"test! thet = "<<thet<<" => "<<gsl_sf_legendre_sphPlm(l,-ml,cos(thet))<<","<<sin(thet)<<std::endl;
      if(thet!=0)
      {
         value=
            (-ml*(pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*cos(-ml*phi))
              +(std::complex<double>(0,-1)*k*nucl_spher_pos[0]*pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*sin(-ml*phi)*dfdf)) // Real spherical harmonics
            *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
      }
      else
      {
         value=
            ((std::complex<double>(0,-1)*k*nucl_spher_pos[0]*pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,-ml,cos(thet))*sin(-ml*phi)*dfdf)) // Real spherical harmonics
            *((pow(std::complex<double>(0,-1),l)*pow(k,l-1)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
      }

   }
   else if(ml>0)
   {
//      std::cout<<"probe 2 f : l = "<<l<<" ; ml = "<<ml<<std::endl;
      if(thet!=0)
      {
         value=
            (-ml*(pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*sin(ml*phi))
            + (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*cos(ml*phi)*dfdf)) // Real spherical harmonics
            *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
            *phase_factor // phase factor
      ;
      }
      else
      {
         value=
            ( (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*pow(-1,ml)*sqrt(2)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*cos(ml*phi)*dfdf)) // Real spherical harmonics
            *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
            *phase_factor // phase factor
            ;
      }
   }
   else
   {
//      std::cout<<"probe 3 f : l = "<<l<<" ; ml = "<<ml<<std::endl;
   value=
      (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*pow(-1,ml)*gsl_sf_legendre_sphPlm(l,ml,cos(thet))*dfdf) // Real spherical harmonics
      *((pow(std::complex<double>(0,-1),l)*pow(k,l)*exp(-pow(k,2)/(4*contraction_zeta)))/(pow(2*contraction_zeta,1.5+l))) // Radial part
      *phase_factor // phase factor
      ;
   }
* /
   return std::conj(value);
}
*/
double MO_value( int mo_index, double x, double y, double z,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size,int **angular_mom_numbers)
{
   double value(0);
   double xp(0);
   double yp(0);
   double zp(0);
   double r(0);
   double thet(0);
   double phi(0);

   for(int i=0;i!=basis_size;i++)
   {
      xp=x-nucl_spher_pos[nucl_basis_func[i]-1][0]*sin(nucl_spher_pos[nucl_basis_func[i]-1][1])*cos(nucl_spher_pos[nucl_basis_func[i]-1][2]);
      yp=y-nucl_spher_pos[nucl_basis_func[i]-1][0]*sin(nucl_spher_pos[nucl_basis_func[i]-1][1])*sin(nucl_spher_pos[nucl_basis_func[i]-1][2]);
      zp=z-nucl_spher_pos[nucl_basis_func[i]-1][0]*cos(nucl_spher_pos[nucl_basis_func[i]-1][1]);
      r=sqrt(pow(xp,2)+pow(yp,2)+pow(zp,2));
      if(r==0)
      {
          thet=acos(-1)/2;
          phi=0;
      }
      else
      {
          thet=acos(zp/r);
      if( xp == 0 && yp >= 0 )
      {
         phi=acos(-1)/2;
      }
      else if(xp == 0 && yp < 0)
      {
         phi=3.*acos(-1)/2; 
      }
      else if(xp < 0 && yp == 0)
      {
         phi=acos(-1); 
      }
      else
      {
         phi=atan2(yp,xp);
      }
   }
      value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_value(i,r,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,basis_func_type,angular_mom_numbers);
//      std::cout<<"==>"<<MO_neut_basis_coeff[mo_index*basis_size+i]<<"*"<<AO_value(i,r,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,basis_func_type,angular_mom_numbers)<<std::endl;
   }
   return value;
}
double AO_value(int ao_index,double r, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int **angular_mom_numbers)
{

   double value(0);
   double norm(0);
   for(int i=0;i!=contraction_number[ao_index];i++)
   {
         value+=contraction_coeff[ao_index][i]*contraction_value(r,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index],angular_mom_numbers[ao_index]);
//         if(value>=1)
//            std::cout<<"==>==>"<<contraction_coeff[ao_index][i]<<" * "<<contraction_value(r,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index],angular_mom_numbers[ao_index])<<" = "<<contraction_coeff[ao_index][i]*contraction_value(r,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index],angular_mom_numbers[ao_index])<<std::endl;
   }
   return value;

}
double contraction_value( double r, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type,int* angular_mom_numbers)
{
   using namespace std;

//   if(exp(-contraction_zeta*pow(r,2))*pow(r,2) > 1)
//      std::cout<<"****"<<exp(-contraction_zeta*pow(r,2))<<","<<r<<std::endl<<"==>"<<exp(-contraction_zeta*pow(r,2))*pow(r,2)<<std::endl;

//   std::cout<<basis_func_type.c_str()<<" - l = "<<angular_mom_numbers[0]<<" ; ml = "<<angular_mom_numbers[1]<<std::endl;

   switch(angular_mom_numbers[0])
   {
      case 0:
         return 0.5*(1/sqrt(acos(-1)))*exp(-contraction_zeta*pow(r,2));
      case 1:
         switch(angular_mom_numbers[1])
         {
            case -1:
               return sqrt(3./(4.*acos(-1)))*sin(thet)*sin(phi)*exp(-contraction_zeta*pow(r,2))*r;
            case 0:
               return sqrt(3./(4.*acos(-1)))*cos(thet)*exp(-contraction_zeta*pow(r,2))*r;
            case 1:
               return sqrt(3./(4.*acos(-1)))*sin(thet)*cos(phi)*exp(-contraction_zeta*pow(r,2))*r;
         }
      case 2:
         switch(angular_mom_numbers[1])
         {
            case -2:
               return 0.25*sqrt(15./(acos(-1)))*pow(sin(thet),2)*sin(2.*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
            case -1:
               return 0.25*sqrt(15./(acos(-1)))*sin(2*thet)*sin(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
            case 0:
               return 0.25*sqrt(5./acos(-1))*(3*pow(cos(thet),2)-1)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
            case 1:
               return 0.25*sqrt(15./(acos(-1)))*sin(2.*thet)*cos(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
            case 2:
               return 0.25*sqrt(15./(acos(-1)))*cos(2*phi)*pow(sin(thet),2)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
         }
      case 3:
         switch(angular_mom_numbers[1])
         {
            case -3:
                return 0.25*sqrt(35/(2*acos(-1)))*(4*pow(cos(phi),2)-1)*pow(sin(thet),3)*sin(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case -2:
                return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*sin(2.*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case -1:
                return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case 0:
                return 0.25*sqrt(7./acos(-1))*cos(thet)*(2*pow(cos(thet),2)-3*pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case 1:
                return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case 2:
                return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*cos(2*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
            case 3:
                return 0.25*sqrt(35/(2*acos(-1)))*(1-4.*pow(sin(phi),2))*pow(sin(thet),3)*cos(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
         }
      default:
         std::cout<<"ERROR : NOT RECOGNIZED ML QUANTUM NUMBER IN CONTRACTION VALUE FUNCTION"<<std::endl;
         exit(EXIT_FAILURE);

   }
   exit(EXIT_FAILURE);
   return 0;
/*
   double temp(0);
   if(basis_func_type=="1s")
   {
      return 0.5*(1/sqrt(acos(-1)))*exp(-contraction_zeta*pow(r,2));
   }
   else if(basis_func_type=="2px")
   {
      return sqrt(3./(4.*acos(-1)))*sin(thet)*cos(phi)*exp(-contraction_zeta*pow(r,2))*r;
   }
   else if(basis_func_type=="2py")
   {
      return sqrt(3./(4.*acos(-1)))*sin(thet)*sin(phi)*exp(-contraction_zeta*pow(r,2))*r;
   }
   else if(basis_func_type=="2pz")
   {
      return sqrt(3./(4.*acos(-1)))*cos(thet)*exp(-contraction_zeta*pow(r,2))*r;
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.25*sqrt(15./(acos(-1)))*pow(sin(thet),2)*sin(2.*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.25*sqrt(15./(acos(-1)))*sin(2*thet)*sin(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
   }
   else if(basis_func_type=="3d0")
   {
      return 0.25*sqrt(5./acos(-1))*(3*pow(cos(thet),2)-1)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.25*sqrt(15./(acos(-1)))*sin(2.*thet)*cos(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.25*sqrt(15./(acos(-1)))*cos(2*phi)*pow(sin(thet),2)*exp(-contraction_zeta*pow(r,2))*pow(r,2);
   }
   else if(basis_func_type=="4f3-")
   {
//      return 0.25*sqrt(35/(2*acos(-1)))*(3*pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*sin(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
      return 0.25*sqrt(35/(2*acos(-1)))*(4*pow(cos(phi),2)-1)*pow(sin(thet),3)*sin(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*sin(2.*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f0")
   {
      return 0.25*sqrt(7./acos(-1))*cos(thet)*(2*pow(cos(thet),2)-3*pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*cos(2*phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.25*sqrt(35/(2*acos(-1)))*(1-4.*pow(sin(phi),2))*pow(sin(thet),3)*cos(phi)*exp(-contraction_zeta*pow(r,2))*pow(r,3);
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
*/
}
bool build_ao_s(double* S,int *nucl_basis_func,int *contraction_number,double **nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,int basis_size,double* lnfact_memo)
{
   double val=0;
   int l_val=0;
   int l_valp=0;
   
   for(int i=0;i!=basis_size;i++)
   {
      l_val=l_number(basis_func_type[i].c_str());
      for(int j=0;j!=basis_size;j++)
      {
         l_valp=l_number(basis_func_type[j].c_str());
         if(nucl_basis_func[i]==nucl_basis_func[j])
         {
            std::cout<<"atom "<<nucl_basis_func[i]<<" - atom "<<nucl_basis_func[j]<<std::endl;
            val=0;
            for(int ip=0;ip!=contraction_number[i];ip++)
            {
               for(int jp=0;jp!=contraction_number[j];jp++)
               {
                  val+=0.5*(contraction_coeff[i][ip]*contraction_coeff[j][jp]*intplushalf_gamma(0.5*l_val+1.5))/(pow(contraction_zeta[i][ip]+contraction_zeta[j][jp],1.5+l_val))*kronecker_delta(l_val,l_valp);
               }
            }
            std::cout<<"<"<<i+1<<"|"<<j+1<<"> = "<<val<<std::endl;
         }
      }
   }
   exit(EXIT_SUCCESS);
   return 1;
}
void build_transition_density_matrix(int n_states_neut,int n_closed,int n_occ,int ci_size_neut,int n_elec_neut,double **ci_vec_neut,double **tran_den_mat_mo)
{

   bool test(0);
   bool test2(0);
   bool test3(0);
   int q(0);
   int p(0);
   double sum(0);
   double det_val;
   
    for (int i=0; i!=n_states_neut; i++)//ELECTRONIC STATE N
    {
        for (int j=0; j!=n_states_neut; j++)//ELECTRONIC STATE K
        {
           std::cout<<" density between states "<<i<<" and "<<j<<std::endl;
           sum=0;
         for(int k=0;k!=(n_closed+n_occ);k++)
         {
            for(int kp=0;kp!=n_closed+n_occ;kp++)
            {
               tran_den_mat_mo[n_states_neut*i+j][(n_occ+n_closed)*k+kp] = 0; 
               for(int m=0;m!=ci_size_neut;m++)
               {
                  for(int n=0;n!=ci_size_neut;n++)
                  {

                     det_val=build_reduced_determinant(k,kp,n_elec_neut,n_closed,n_occ,&ci_vec_neut[0][(n_elec_neut+n_states_neut)*m],&ci_vec_neut[0][(n_elec_neut+n_states_neut)*n],&ci_vec_neut[1][n_elec_neut*m],&ci_vec_neut[1][n_elec_neut*n]);

                     tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]+=ci_vec_neut[0][(n_elec_neut+n_states_neut)*(m)+n_elec_neut+i]*ci_vec_neut[0][(n_elec_neut+n_states_neut)*(n)+n_elec_neut+j]*det_val;
                    }
                 }
               if(k==kp)
               {
                 sum+=tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp];
               //  std::cout<<" from orbital "<<k<<" and from orbital "<<kp<<":"<<tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]<<std::endl;
               }
//               std::cout<<std::setprecision(8)<<"trdm val "<<tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]<<std::endl;
              }
           }
         std::cout<<"SUM = "<<sum<<std::endl;
        }
    }
}
double build_reduced_determinant( int ai,int aj,int n_elec,int n_closed,int n_occ,double* mo_vector_1,double* mo_vector_2,double *spin_vector_1,double *spin_vector_2)
{
   /* Given the vectors containing the mo labels and the spin labels of the electrons, this routine builds a slater determinant from which one electron contained in the mo's i and j have been removed   !!!! ONLY FOR SINGLET AND SIMPLE EXCITATION
   */

   bool test2(0);
   bool test3(0);
   int spin(0);
   double temp(0);

   int new_vector_1[(n_occ+n_closed)];
   int new_vector_2[(n_occ+n_closed)];

   double prefactor(1);

   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      new_vector_1[k]=0;
      new_vector_2[k]=0;
   }
   
   for(int e=0;e!=n_elec;e++)
   {
      new_vector_1[int(mo_vector_1[e])]+=1;
      new_vector_2[int(mo_vector_2[e])]+=1;
   }
   /*
   std::cout<<"Taking electron from orbitals "<<ai<<","<<aj<<std::endl;
   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      std::cout<<new_vector_1[k]<<" ";
   }std::cout<<std::endl;
   
   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      std::cout<<new_vector_2[k]<<" ";
   }std::cout<<std::endl;
   */
   prefactor=sqrt(double(new_vector_1[ai]))*sqrt(double(new_vector_2[aj]));
   new_vector_1[ai]--;
   new_vector_2[aj]--;

   for(int k=0;k!=(n_occ+n_closed);k++)
   {
//      std::cout<<new_vector_1[k]<<std::endl;
      if(new_vector_1[k]!=new_vector_2[k])
         return 0;
   }
   
  // if(prefactor != 0 )
  //    std::cout<<prefactor<<std::endl;
   return prefactor;
}
void compute_bessel_pice_mo(double*** pice_ortho_mo,double*** pice_ddx_mo,double*** pice_ddy_mo,double*** pice_ddz_mo,int jl_max,int n_occ,int basis_size,int nk,double kmax,double *MO_coeff_neutral,double **contraction_zeta,double **contraction_coeff,int * contraction_number,double** nucl_spher_pos,int *nucl_basis_func,int** angular_mom_numbers)
{

   int ll2(0);
   int mm2(0);

   int l2max(0);
   int max_contraction_num(0);

   double kp(0);

   for(int ww=0;ww!=basis_size;ww++)
   {
      if(l2max<angular_mom_numbers[ww][0])
         l2max=angular_mom_numbers[ww][0];
      if(max_contraction_num<contraction_number[ww])
         max_contraction_num=contraction_number[ww];
   }
   int l3max=jl_max+l2max+1;

   std::cout<<"ALLOCATING BESSEL DERIVED ARRAYS"<<std::endl;
   double ***Kval_cont=new double **[basis_size];
   double ***ddk_Kval_cont=new double **[basis_size];
   double **Kval=new double *[basis_size];
   double **ddk_Kval=new double *[basis_size];
   double ***bessel_val=new double **[basis_size];
   double ***ddk_bessel_val=new double **[basis_size];
   //There are 9 angular integrals to compute.
   //for a given value of l_max, there are in total S_{l_max}=∑_l=0^{l_max} 2l+1 values of ml
   //S_{l_max}= 2 * 0.5 * l_max * (l_max + 1) + l_max + 1 = l_max^2 + 2 * l_max + 1
   double ***ang_int1=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int2=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int3=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int4=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int5=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int6=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int7=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int8=new double ** [jl_max*jl_max+2*jl_max+1]; 
   double ***ang_int9=new double ** [jl_max*jl_max+2*jl_max+1]; 
   //Then, we need the different cartesian components of the pice

   double ***pice_ddx_basis=new double **[jl_max*jl_max+2*jl_max+1];
   double ***pice_ddy_basis=new double **[jl_max*jl_max+2*jl_max+1];
   double ***pice_ddz_basis=new double **[jl_max*jl_max+2*jl_max+1];
   double ***pice_ortho_basis=new double **[jl_max*jl_max+2*jl_max+1];

   for( int ww = 0 ; ww != basis_size; ww ++)
   {
      Kval_cont[ww]=new double *[max_contraction_num];
      ddk_Kval_cont[ww]=new double *[max_contraction_num];

      Kval[ww]=new double [nk];
      ddk_Kval[ww]=new double [nk];

      bessel_val[ww]=new double *[l3max*l3max+2*l3max+1];
      ddk_bessel_val[ww]=new double *[l3max*l3max+2*l3max+1];

      for(int bb=0;bb!=max_contraction_num;bb++)
      {
         Kval_cont[ww][bb]=new double [nk];
         ddk_Kval_cont[ww][bb]=new double [nk];
      }
      for(int bb=0;bb!=l3max*l3max+2*l3max+1;bb++)
      {
         bessel_val[ww][bb]=new double [nk];
         ddk_bessel_val[ww][bb]=new double [nk];
      }
   }
   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   {
      ang_int1[ji]=new double * [basis_size];
      ang_int2[ji]=new double * [basis_size];
      ang_int3[ji]=new double * [basis_size];
      ang_int4[ji]=new double * [basis_size];
      ang_int5[ji]=new double * [basis_size];
      ang_int6[ji]=new double * [basis_size];
      ang_int7[ji]=new double * [basis_size];
      ang_int8[ji]=new double * [basis_size];
      ang_int9[ji]=new double * [basis_size];
      pice_ddx_basis[ji]=new double *[basis_size];
      pice_ddy_basis[ji]=new double *[basis_size];
      pice_ddz_basis[ji]=new double *[basis_size];
      pice_ortho_basis[ji]=new double *[basis_size];
      for( int ww = 0 ; ww != basis_size; ww ++)
      {
         ang_int1[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int2[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int3[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int4[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int5[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int6[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int7[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int8[ji][ww]=new double[l3max*l3max+2*l3max+1];
         ang_int9[ji][ww]=new double[l3max*l3max+2*l3max+1];
         pice_ddx_basis[ji][ww]=new double [nk];
         pice_ddy_basis[ji][ww]=new double [nk];
         pice_ddz_basis[ji][ww]=new double [nk];
         pice_ortho_basis[ji][ww]=new double [nk];
      }
   }
   std::cout<<"BESSEL DERIVED ARRAYS ALLOCATED!!"<<std::endl;
   for(int ww=0;ww!=basis_size;ww++)//l2
   {
      ll2=angular_mom_numbers[ww][0];
      mm2=angular_mom_numbers[ww][1];
      std::cout<<"Basis function "<<ww<<"/"<<basis_size<<std::endl;

      for(int ll3=0;ll3!=l3max+1;ll3++)//l3
      {
         for(int mm3=-ll3;mm3!=ll3+1;mm3++)
         {
            for(int ll1=0;ll1!=jl_max+1;ll1++)//l1
            {
               for(int mm1=-ll1;mm1!=ll1+1;mm1++)
               {
                      ang_int1[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);
                      ang_int2[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);

                      ang_int3[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_D_integral(mm1,mm2,mm3);

                      ang_int4[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

                      ang_int5[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

                      ang_int6[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_D_integral(mm1,mm2,mm3);

                      ang_int7[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);

                      ang_int8[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);
                      ang_int9[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]][1],nucl_spher_pos[nucl_basis_func[ww]][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *gaunt_formula(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);
               }
            }
         }
       }

       for(int k=0;k!=nk;k++)
       {
            //Computing the k-dependent part of the integral
            //there are two types of integral to compute. the direct expresison and the derivative of the expression
            //This part only depends on l2 and l3
         kp=kmax*(k+1)/nk;
         Kval[ww][k]=0;
         ddk_Kval[ww][k]=0;

         for(int bb=0;bb!=contraction_number[ww];bb++)
         {
               Kval_cont[ww][bb][k]=pow(kp,ll2)*exp(-kp*kp/(4*contraction_zeta[ww][bb]))/pow(2*contraction_zeta[ww][bb],1.5+ll2);
               ddk_Kval_cont[ww][bb][k]=(double(ll2)/kp-kp/(2*contraction_zeta[ww][bb]))*Kval_cont[ww][bb][k];

               Kval[ww][k]+=contraction_coeff[ww][bb]*Kval_cont[ww][bb][k];
               ddk_Kval[ww][k]+=contraction_coeff[ww][bb]*ddk_Kval_cont[ww][bb][k];
         }

         for(int ll3=0;ll3!=l3max+1;ll3++)//l3
         {
            bessel_val[ww][ll3][k]=j_l(ll3,kp*nucl_spher_pos[nucl_basis_func[ww]][0]);
            ddk_bessel_val[ww][ll3][k]=nucl_spher_pos[nucl_basis_func[ww]][0]*dj_ldz(ll3,kp*nucl_spher_pos[nucl_basis_func[ww]][0]);
         }
       }
   }

//Now we have computed all the angular and radial integrals for the basis functions. We have as much angular integrals computed as basis functions, Bessel functions and plane wave spherical expansion functions.
//We need to combine these integrals in order to get the values for the MO. 
//First, for a single l1,m2,l2,m2, we need to sum up the functions for the values of l3,m3
//
  int l1(0);
  bool test(0);

   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   {
      test=0;
      while(!test)
      {
         if(ji<=l1*l1+2*l1)
            break;
         l1++;
      }
      
      for( int ww = 0 ; ww != basis_size; ww ++)
      {
        ll2=angular_mom_numbers[ww][0];
        mm2=angular_mom_numbers[ww][1];
         for( int k=0;k!=nk;k++)
         {
            kp=kmax*(k+1)/nk;
            pice_ddx_basis[ji][ww][k]=0;
            pice_ddy_basis[ji][ww][k]=0;
            pice_ddz_basis[ji][ww][k]=0;
            pice_ortho_basis[ji][ww][k]=0;
            for(int ll3=0;ll3!=l3max+1;ll3++)//l3
            {
               for(int mm3=-ll3;mm3!=ll3+1;mm3++)
               {
                  pice_ddx_basis[ji][ww][k]+=
                     (kp*(ddk_Kval[ww][k]*bessel_val[ww][ll3][k]+Kval[ww][k]*ddk_bessel_val[ww][ll3][k])*ang_int1[ji][ww][ll3*ll3+ll3+mm3]
                     +Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int2[ji][ww][ll3*ll3+ll3+mm3]
                     +Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int3[ji][ww][ll3*ll3+ll3+mm3]);

                  pice_ddy_basis[ji][ww][k]+=
                     (kp*(ddk_Kval[ww][k]*bessel_val[ww][ll3][k]+Kval[ww][k]*ddk_bessel_val[ww][ll3][k])*ang_int4[ji][ww][ll3*ll3+ll3+mm3]
                     +Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int5[ji][ww][ll3*ll3+ll3+mm3]
                     +Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int6[ji][ww][ll3*ll3+ll3+mm3]);

                  pice_ddz_basis[ji][ww][k]+=
                     (kp*(ddk_Kval[ww][k]*bessel_val[ww][ll3][k]+Kval[ww][k]*ddk_bessel_val[ww][ll3][k])*ang_int7[ji][ww][ll3*ll3+ll3+mm3]
                     +Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int8[ji][ww][ll3*ll3+ll3+mm3]);

                  pice_ortho_basis[ji][ww][k]+=
                     (kp*Kval[ww][k]*bessel_val[ww][ll3][k]*ang_int9[ji][ww][ll3*ll3+ll3+mm3]);
               }
            }
         }
      }
   }
   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   {
      for(int mm=0;mm!=n_occ;mm++)
      {
         for( int k=0;k!=nk;k++)
         {
            pice_ortho_mo[ji][mm][k]=0;
            pice_ddx_mo[ji][mm][k]=0;
            pice_ddy_mo[ji][mm][k]=0;
            pice_ddz_mo[ji][mm][k]=0;
            for( int ww = 0 ; ww != basis_size; ww ++)
            {
               pice_ddx_mo[ji][mm][k]+=MO_coeff_neutral[mm*basis_size+ww]*pice_ddx_basis[ji][ww][k];
               pice_ddy_mo[ji][mm][k]+=MO_coeff_neutral[mm*basis_size+ww]*pice_ddy_basis[ji][ww][k];
               pice_ddz_mo[ji][mm][k]+=MO_coeff_neutral[mm*basis_size+ww]*pice_ddz_basis[ji][ww][k];
               pice_ortho_mo[ji][mm][k]+=MO_coeff_neutral[mm*basis_size+ww]*pice_ortho_basis[ji][ww][k];
            }
         }
      }
   }

   delete [] pice_ddx_basis;
   delete [] pice_ddy_basis;
   delete [] pice_ddz_basis;
   delete [] pice_ortho_basis;

   delete [] Kval_cont; 
   delete [] ddk_Kval_cont; 
   delete [] Kval; 
   delete [] ddk_Kval; 
   delete [] bessel_val;
   delete [] ddk_bessel_val;
   delete [] ang_int1;
   delete [] ang_int2;
   delete [] ang_int3;
   delete [] ang_int4;
   delete [] ang_int5;
   delete [] ang_int6;
   delete [] ang_int7;
   delete [] ang_int8;
   delete [] ang_int9;
}
