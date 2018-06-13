bool dyson_mo_coeff_comp(int n_states_neut,int n_states_cat,int n_occ,int ci_size_neut,int ci_size_cat,int n_elec_neut,double **ci_vec_neut,double **ci_vec_cat,double *overlap,double *Dyson_MO_basis_coeff)
{
   bool test(0);
   bool test2(0);

   double *temp;
   temp=new double[(n_elec_neut)*(n_elec_neut)];

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
                    //std::cout<<"New configuration of the neutral : "<<n<<std::endl;
                    for (int l=0; l!=ci_size_cat; l++)//  over configuration of the cation
                    {
                        test2=0;
                     //========================VVVVVVVVVVVVV Overlap matrix for a given configuration VVVVVVVVVVVVVVVVV==================   
                        for (int m=0; m!=n_elec_neut; m++)  //Over the electrons of the neutral
                        {
                            for (int p=0; p!=n_occ; p++)//Over the MO of the neutral
                            {
                                    if(int(ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+m])==p)
                                    {
                                        if (p==k && ci_vec_neut[1][(n_elec_neut)*n+m])
                                        {
                                            //std::cout<<" p = "<<p<<" k = "<<k<<" =>  taking electron ß"<<std::endl;
                                            test=1;
                                            test2=1;
                                            continue;
                                        }
                                        
                                        for (int o=0; o!=n_elec_neut-1; o++)//Over the electrons of the cation
                                        {
                                            for (int q=0; q!=n_occ; q++)//Over the MO of the cation
                                            {
                                                if(int(ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+o])==q)
                                                {
                                                    temp[(n_elec_neut-1)*(m-test2)+o]=overlap[n_occ*p+q]*kronecker_delta(ci_vec_neut[1][n_elec_neut*n+(m-test2)], ci_vec_cat[1][(n_elec_neut-1)*l+o]);
                                                }
                                            }
                                        }
                                        
                                    }
                            }
                        }
                        //========================^^^^^^^^^^^^^^ Overlap matrix for a given configuration ^^^^^^^^^^^^^==================
                        
                        
/*                        //====================CHECKING DETERMINANT===============
                        
                        for (int m=0; m!=n_elec_neut; m++)
                        {
                            std::cout<<(ci_vec_neut[1][n_elec_neut*n+m])<<"    ";
                            
                        }std::cout<<std::endl;
                        for (int m=0; m!=n_elec_neut-1; m++)
                        {
                            std::cout<<(ci_vec_cat[1][(n_elec_neut-1)*l+m])<<"    ";
                            
                        }std::cout<<std::endl;
                        for (int m=0; m!=n_elec_neut-1; m++)
                        {
                            for (int o=0; o!=n_elec_neut-1; o++)
                            {
                                std::cout<<std::setw(12)<<std::setprecision(5)<<temp[(n_elec_neut-1)*m+o];
                            }std::cout<<std::endl;
                        }std::cout<<std::endl;
                        
                        
                        
                        //====================CHECKING DETERMINANT===============*/
                        if(test2)
                        {
                        Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]+=(ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]*ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]*determinant(temp,(n_elec_neut-1)));
                           //std::cout<<" sum is "<<Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]<<std::endl;
                        }

                        
                        //std::cout<<std::endl<<std::setw(12)<<std::setprecision(5)<<"config neutral "<<n<<" = "<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]<<std::setw(12)<<std::setprecision(5)<<"  config cation "<<l<<" = "<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]<<std::setw(12)<<std::setprecision(5)<<determinant(temp,(n_elec_neut-1))<<std::endl<<std::endl;
                        
                        //std::cout<<"config neutral "<<n<<" = "<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]<<"  config cation "<<l<<" = "<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]<<std::endl;
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

std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type)
{
   return 0;
}
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type)
{
   using namespace std;

   complex<double> phase_factor(exp(-std::complex<double>(0,k*nucl_spher_pos[0]*(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1])))));
   if(basis_func_type=="1s")
   {
      /*
       * 1/2 Sqrt[1/Pi] (1/(2Pi))^(3/2) Exp[-I k r f[\[Theta],\[Phi],\[CapitalTheta],\[CapitalPhi]]]E^(-(k^2/(4 d)))
       * */
      return 0.017911224007836134*exp(-k*k/(4*contraction_zeta))*phase_factor;
   }
   else if(basis_func_type=="2px")
   {
      /*
       * Sqrt[3/(4Pi)]Sin[\[Theta]]Sin[\[Phi]]4Pi I Exp[-I k r  f[\[Theta],\[Phi],\[CapitalTheta],\[CapitalPhi]]] ((-(1/d))^(3/2) k (-Sqrt[d] k+(2 d+k^2) DawsonF[k/(2 Sqrt[d])]))/(2 (-(k^2/d))^(3/2) \[Pi]^(3/2))
       * */
      return 0.5513288954217921*sin(thet)*cos(phi)*phase_factor*(-sqrt(contraction_zeta)+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))));
   }
   else if(basis_func_type=="2py")
   {
      return 0.5513288954217921*sin(thet)*sin(phi)*phase_factor*(-sqrt(contraction_zeta)+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))));
   }
   else if(basis_func_type=="2pz")
   {
      return 0.5513288954217921*cos(thet)*phase_factor*(-sqrt(contraction_zeta)+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))));
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.6307831305050401*phase_factor*sin(thet)*sin(thet)*sin(phi)*cos(phi)*(-0.6752372371178296*pow(contraction_zeta,1.5)*erf(k/(2*contraction_zeta))+k*exp(-k*k/4*contraction_zeta)*(6*contraction_zeta+k*k))/pow(k,3);
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.6307831305050401*phase_factor*sin(thet)*cos(theta)*sin(phi)*(-0.6752372371178296*pow(contraction_zeta,1.5)*erf(k/(2*contraction_zeta))+k*exp(-k*k/4*contraction_zeta)*(6*contraction_zeta+k*k))/pow(k,3);
   }
   else if(basis_func_type=="3d0")
   {
      return 0.04005071444174783*phase_factor*(-sin(thet)*sin(thet)*cos(phi)*cos(phi)-sin(thet)*sin(thet)*sin(phi)*sin(phi)+2*cos(theta)*cos(theta))*(-0.6752372371178296*pow(contraction_zeta,1.5)*erf(k/(2*contraction_zeta))+k*exp(-k*k/4*contraction_zeta)*(6*contraction_zeta+k*k))/pow(k,3);
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.6307831305050401*phase_factor*sin(thet)*cos(thet)*cos(phi)*(-0.6752372371178296*pow(contraction_zeta,1.5)*erf(k/(2*contraction_zeta))+k*exp(-k*k/4*contraction_zeta)*(6*contraction_zeta+k*k))/pow(k,3);
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.3376186185589148*phase_factor*(sin(thet)*sin(thet)*cos(phi)*cos(phi)-sin(thet)*sin(thet)*sin(phi)*sin(phi))*(-0.6752372371178296*pow(contraction_zeta,1.5)*erf(k/(2*contraction_zeta))+k*exp(-k*k/4*contraction_zeta)*(6*contraction_zeta+k*k))/pow(k,3);
   }
   else if(basis_func_type=="4f3-")
   {}
   else if(basis_func_type=="4f2-")
   {}
   else if(basis_func_type=="4f1-")
   {}
   else if(basis_func_type=="4f0")
   {}
   else if(basis_func_type=="4f1+")
   {}
   else if(basis_func_type=="4f2+")
   {}
   else if(basis_func_type=="4f3+")
   {}
   else
   {
      exit(EXIT_FAILURE);
   }
}
