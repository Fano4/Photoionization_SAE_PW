void compute_bessel_pice_mo(double*** pice_ortho_mo,double*** pice_ddx_mo,double*** pice_ddy_mo,double*** pice_ddz_mo,int jl_max,int n_occ,int basis_size,int nk,double kmax,double *MO_coeff_neutral,double **contraction_zeta,double **contraction_coeff,int * contraction_number,double** nucl_spher_pos,int *nucl_basis_func,int** angular_mom_numbers)
{

   int ll2(0);
   int mm2(0);
   int l2max(0);
   int max_contraction_num(0);
   double antiphase(1); //THE GAUNT INTEGRAL AND RECURRENCE RELATIONS ARE BASED ON EXPRESSIONS INCLUDING THE CONDON-SHORTLEY PHASE
                        // WHICH IS EXPLICITELY ACCOUNTED FOR IN THE SPHERICAL HARMINCS PREFACTORS. TO AVOIS COUNTING IT TWICE, 
                        // WE ADD AN ANTIPHASE IN THE INTEGRAL EXPRESSIONS
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
                  antiphase=pow(-1,mm1+mm2+mm3);
//                  std::cout<<"l1 = "<<ll1<<", l2 = "<<ll2<<", l3 = "<<ll3<<", m1 = "<<mm1<<", m2 = "<<mm2<<", m3 = "<<mm3<<std::endl;
//                  std::cout<<"ang_int1 , ";
                      ang_int1[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int2 , ";
                      ang_int2[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int3 , ";
                      ang_int3[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_D_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int4 , ";
                      ang_int4[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int5 , ";
                      ang_int5[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_m1_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int6 , ";
                      ang_int6[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m2(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*I_p1_D_integral(mm1,mm2,mm3);

//                  std::cout<<"ang_int7 , ";
                      ang_int7[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_p1(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);

//                  std::cout<<"ang_int8 , ";
                      ang_int8[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1-1)/2))
                          *(-4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *J_int_m1_D(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);

//                  std::cout<<"ang_int9 , ";
                      ang_int9[ll1*ll1+ll1+mm1][ww][ll3*ll3+ll3+mm3]=
                         antiphase*pow(-1,((ll2+ll3-ll1)/2))
                          *(4.*acos(-1)*rYlm(ll3,mm3,nucl_spher_pos[nucl_basis_func[ww]-1][1],nucl_spher_pos[nucl_basis_func[ww]-1][2]))
                          *prefactor_rYlm(ll1,fabs(mm1))*prefactor_rYlm(ll2,fabs(mm2))*prefactor_rYlm(ll3,fabs(mm3))
                          *gaunt_formula(ll1,ll2,ll3,fabs(mm1),fabs(mm2),fabs(mm3))*azim_integ(mm1,mm2,mm3);
//                      std::cout<<"END"<<std::endl;
//                      exit(EXIT_SUCCESS);
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
