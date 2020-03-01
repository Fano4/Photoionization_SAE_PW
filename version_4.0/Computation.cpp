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
//                  std::cout<<"l1 = "<<ll1<<", l2 = "<<ll2<<", l3 = "<<ll3<<", m1 = "<<mm1<<", m2 = "<<mm2<<", m3 = "<<mm3<<std::endl;
//                  std::cout<<"ang_int1 , ";
                      ang_int1[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);

//                  exit(EXIT_SUCCESS);
//                  std::cout<<"ang_int2 , ";
                      ang_int2[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int3 , ";
                      ang_int3[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_D_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int4 , ";
                      ang_int4[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int5 , ";
                      ang_int5[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int6 , ";
                      ang_int6[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_D_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int7 , ";
                      ang_int7[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);

//                  std::cout<<"ang_int8 , ";
                      ang_int8[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);

//                  std::cout<<"ang_int9 , ";
                      ang_int9[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         pow(-1,((ll2+ll3-ll1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *gaunt_formula(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);
//                      std::cout<<"END"<<std::endl;
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
            bessel_val[ww][ll3][k]=j_l(ll3,kp*nucl_spher_pos[nucl_basis_func[ww]-1][0]);
            ddk_bessel_val[ww][ll3][k]=nucl_spher_pos[nucl_basis_func[ww]-1][0]*dj_ldz(ll3,kp*nucl_spher_pos[nucl_basis_func[ww]-1][0]);
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
