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

std::complex<double> MO_Fourier_transform_grad( int mo_index,int comp, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double **MO_neut_basis_coeff,int basis_size)
{
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT_grad(i,comp,k,thet,phi,contraction_number,nucl_spher_pos,nucl_basis_func,contraction_coeff,contraction_zeta,basis_func_type);
   }
   return value;
}
std::complex<double> AO_FT_grad(int ao_index,int comp,double k, double thet, double phi,int *contraction_number,double **nucl_spher_pos,int *nucl_basis_func,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type)
{
   std::complex<double> value(0);
   switch(comp)
   {
      case 0:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FTi_grad_k(k,thet,phi,nucl_sphere_pos[ao_index],contraction_zeta[ao_index][i],nucl_basis_func[ao_index]);
       }
       return value;
      case 1:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FT_grad_thet(k,thet,phi,nucl_sphere_pos[ao_index],contraction_zeta[ao_index][i],nucl_basis_func[ao_index]);
       }
       return value;
      case 2:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FT_grad_phi(k,thet,phi,nucl_sphere_pos[ao_index],contraction_zeta[ao_index][i],nucl_basis_func[ao_index]);
       }
       return value;
      default:
       std::cout<<"invalid component of grad vector required"<<std::endl;
       exit(EXIT_FAILURE);
   }

}
std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double **MO_neut_basis_coeff,int basis_size)
{
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT(i,k,thet,phi,contraction_number,nucl_spher_pos,nucl_basis_func,contraction_coeff,contraction_zeta,basis_func_type);
   }
   return value;
}
std::complex<double> AO_FT(int ao_index,double k, double thet, double phi,int *contraction_number,double **nucl_spher_pos,int *nucl_basis_func,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type)
{
   std::complex<double> value(0);
   for(int i=0;i!=contraction_number[ao_index];i++)
   {
      value+=contraction_coeff[ao_index][i]*contraction_FT(k,thet,phi,nucl_sphere_pos[ao_index],contraction_zeta[ao_index][i],nucl_basis_func[ao_index]);
   }
   return value;

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
   {
      return std::complex<double(0,1)*0.042273611654255666*(3*pow(sin(thet)*cos(phi),2)-pow(sin(thet)*sin(phi),2))*sin(thet)*sin(phi)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f2-")
   {
      return std::complex<double(0,1)*0.20709755627499729*(pow(sin(thet),2)*cos(thet)*sin(phi)*cos(phi)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f1-")
   {
      return std::complex<double(0,1)*0.046308421380498219*(4*pow(cos(thet),2)-pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*sin(phi)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f0")
   {
      return std::complex<double(0,1)*0.026736179549777268*(2*pow(cos(thet),2)-3*pow(sin(thet),2)*pow(cos(phi),2)-3*pow(sin(thet),2)*pow(sin(phi),2))*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f1+")
   {
      return std::complex<double(0,1)*0.046308421380498219*(4*pow(cos(thet),2)-pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*cos(phi)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f2+")
   {
      return std::complex<double(0,1)*0.20709755627499729*cos(thet)*((pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(phi),2)*pow(sin(thet),2)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else if(basis_func_type=="4f3+")
   {
      return std::complex<double(0,1)*0.042273611654255666*(pow(sin(thet)*cos(phi),2)-pow(sin(thet)*sin(phi),2))*sin(thet)*sicos(phi)*phase_factor*((60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k));
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
}
std::complex<double> contraction_FT_grad_k( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type)
{
   using namespace std;

   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1])));
   complex<double> phase_factor(exp(-std::complex<double>(0,k*nucl_spher_pos[0]*f)));
   if(basis_func_type=="1s")
   {
      return 0.0089556120039180672*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k-2*std::complex<double>(0,1)*contraction_zeta*nucl_spher_pos[0]*f)
   }
   else if(basis_func_type=="2px")
   {
      return 0.27566444771089604*phase_factor*sin(thet)*cos(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(4*contraction_zeta+k*k)-2*contraction_zeta*k*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(-std::complex<double>(0,1)*(8*contraction_zeta*contraction_zeta+2*contraction_zeta*k*k+k*k*k*k)+2*contraction_zeta*k*(2*contraction_zeta+k*k)*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k);
   }
   else if(basis_func_type=="2py")
   {
      return 0.27566444771089604*phase_factor*sin(thet)*sin(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(4*contraction_zeta+k*k)-2*contraction_zeta*k*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(-std::complex<double>(0,1)*(8*contraction_zeta*contraction_zeta+2*contraction_zeta*k*k+k*k*k*k)+2*contraction_zeta*k*(2*contraction_zeta+k*k)*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k);
   }
   else if(basis_func_type=="2pz")
   {
      return 0.27566444771089604*phase_factor*cos(thet)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(4*contraction_zeta+k*k)-2*contraction_zeta*k*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(-std::complex<double>(0,1)*(8*contraction_zeta*contraction_zeta+2*contraction_zeta*k*k+k*k*k*k)+2*contraction_zeta*k*(2*contraction_zeta+k*k)*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k);
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.034684936146269905*cos(phi)*sin(thet)*sin(thet)*sin(phi)*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k*(36*contraction_zeta*contraction_zeta+6*contraction_zeta*k*k+k*k*k*k)-2*std::complex<double>(0,1)*contraction_zeta*k*k*(6*contraction_zeta+k*k)*nucl_spher_pos[0]*f+12*pow(contraction_zeta,2.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta)))*(3+std::complex<double>(0,1)*k*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k*k);
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.034684936146269905*cos(thet)*sin(thet)*sin(phi)*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k*(36*contraction_zeta*contraction_zeta+6*contraction_zeta*k*k+k*k*k*k)-2*std::complex<double>(0,1)*contraction_zeta*k*k*(6*contraction_zeta+k*k)*nucl_spher_pos[0]*f+12*pow(contraction_zeta,2.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta)))*(3+std::complex<double>(0,1)*k*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k*k);
   }
   else if(basis_func_type=="3d0")
   {
      return 0.0050063393052184784*(1+3*cos(2*thet))*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k*(36*contraction_zeta*contraction_zeta+6*contraction_zeta*k*k+k*k*k*k)-2*std::complex<double>(0,1)*contraction_zeta*k*k*(6*contraction_zeta+k*k)*nucl_spher_pos[0]*f+12*pow(contraction_zeta,2.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta)))*(3+std::complex<double>(0,1)*k*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k*k);
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.034684936146269905*cos(thet)*sin(thet)*cos(phi)*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k*(36*contraction_zeta*contraction_zeta+6*contraction_zeta*k*k+k*k*k*k)-2*std::complex<double>(0,1)*contraction_zeta*k*k*(6*contraction_zeta+k*k)*nucl_spher_pos[0]*f+12*pow(contraction_zeta,2.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta)))*(3+std::complex<double>(0,1)*k*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k*k);
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.017342468073134953*cos(2*phi)*sin(thet)*sin(thet)*phase_factor*exp(-k*k/(4*contraction_zeta))*(-k*(36*contraction_zeta*contraction_zeta+6*contraction_zeta*k*k+k*k*k*k)-2*std::complex<double>(0,1)*contraction_zeta*k*k*(6*contraction_zeta+k*k)*nucl_spher_pos[0]*f+12*pow(contraction_zeta,2.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta)))*(3+std::complex<double>(0,1)*k*f*nucl_spher_pos[0]))/(contraction_zeta*k*k*k*k);
   }
   else if(basis_func_type=="4f3-")
   {
      return 0.021136805827127833*phase_factor*(1+cos(2*phi))*pow(sin(thet),3)*sin(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.10354877813749866*phase_factor*cos(thet)*cos(phi)*pow(sin(thet),2)*sin(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.0081862496960485968*phase_factor*(3+5*cos(2*thet))*sin(thet)*sin(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f0")
   {
      return 0.006684044887444316*phase_factor*(-1+5*cos(2*thet))*cos(thet)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.0081862496960485968*phase_factor*(3+5*cos(2*phi))*sin(thet)*cos(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.10354877813749866*phase_factor*cos(thet)*cos(2*phi)*pow(sin(thet),2)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.021136805827127833*phase_factor*(-1+2*cos(2*phi))*pow(sin(thet),3)*cos(phi)*(sqrt(contraction_zeta)*k*(std::complex<double>(0,1)*(6*contraction_zeta*k+k*k*k)-2*contraction_zeta*(30*contraction_zeta+k*k)*f*nucl_spher_pos[0])+gsl_sf_dawson(k/(2*sqrt(contraction_zeta)))*(std::complex<double>(0,-1)*(12*contraction_zeta*contraction_zeta*k+4*contraction_zeta*k*k*k+pow(k,5))+2*contraction_zeta*(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*nucl_spher_pos[0]*f))/contraction_zeta;
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
}
std::complex<double> contraction_FT_grad_thet( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type)
{
   using namespace std;

   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdt(cos(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])-sin(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdf(-sin(thet)*sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(-std::complex<double>(0,k*nucl_spher_pos[0]*f)));

   if(basis_func_type=="1s")
   {
      return 0.017911224007836134*std::complex<double>(0,-1)*phase_factor*exp(-k*k/(4*contraction_zeta))*k*nucl_spher_pos[0]*dfdt;
   }
   else if(basis_func_type=="2px")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*cos(phi)*(std::complex<double>(0,1)*cos(thet)+k*nucl_spher_pos[0]*sin(theta)*dfdt)/(k*k);
   }
   else if(basis_func_type=="2py")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(phi)*(std::complex<double>(0,1)*cos(thet)+k*nucl_spher_pos[0]*sin(theta)*dfdt)/(k*k);
   }
   else if(basis_func_type=="2pz")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(std::complex<double>(0,-1)*sin(thet)+k*nucl_spher_pos[0]*cos(theta)*dfdt)/(k*k);
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.069369872292539811*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*cos(phi)*sin(phi)*sin(thet)*(2*cos(thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(thet)*dfdt)/pow(k,3);
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.034684936146269905*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*sin(phi)*(2*cos(2*thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(2*thet)*dfdt)/pow(k,3);
   }
   else if(basis_func_type=="3d0")
   {
      return 0.010012678610436957*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*(6*sin(2*thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*(1+3*cos(2*thet))*dfdt)/pow(k,3);
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.034684936146269905*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*cos(phi)*(2*cos(2*thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(2*thet)*dfdt)/pow(k,3);
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.034684936146269905*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*cos(2*phi)*sin(thet)*(2*cos(thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(thet)*dfdt)/pow(k,3);
   }
   else if(basis_func_type=="4f3-")
   {
      return 0.042273611654255666*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(1+2*cos(2*phi))*sin(thet)*sin(thet)*sin(phi)*(3*std::complex<double>(0,1)*cos(thet)+k*nucl_spher_pos[0]*sin(thet)*dfdt);
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.05177438906874933*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*sin(2*phi)*(std::complex<double>(0,1)+3*std::complex<double>(0,1)*cos(2*thet)+k*nucl_spher_pos[0]*sin(2*thet)*dfdt);
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.0081862496960485968*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(phi)*(std::complex<double>(0,1)*(cos(thet)+15*cos(3*thet))+k*nucl_spher_pos[0]*(sin(thet)*+5*sin(3*thet))*dfdt);
   }
   else if(basis_func_type=="4f0")
   {
      return 0.006684044887444316*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(-3*std::complex<double>(0,1)*(sin(thet)+5*sin(3*thet))+k*nucl_spher_pos[0]*(3*cos(thet)*+5*cos(3*thet))*dfdt);
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.0081862496960485968*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*cos(phi)*(std::complex<double>(0,1)*(cos(thet)+15*cos(3*thet))+k*nucl_spher_pos[0]*(sin(thet)*+5*sin(3*thet))*dfdt);
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.05177438906874933*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*cos(2*phi)*(std::complex<double>(0,1)+3*std::complex<double>(0,1)*cos(2*thet)+k*nucl_spher_pos[0]*sin(2*thet)*dfdt);
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.042273611654255666*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(1+2*cos(2*phi))*sin(thet)*sin(thet)*cos(phi)*(3*std::complex<double>(0,1)*cos(thet)+k*nucl_spher_pos[0]*sin(thet)*dfdt);
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
}
std::complex<double> contraction_FT_grad_phi( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type)
{
   using namespace std;

   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdt(cos(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])-sin(thet)*cos(nucl_spher_pos[1]));
   complex<double> dfdf(-sin(thet)*sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(-std::complex<double>(0,k*nucl_spher_pos[0]*f)));
   if(basis_func_type=="1s")
   {
      return 0.017911224007836134*std::complex<double>(0,-1)*phase_factor*exp(-k*k/(4*contraction_zeta))*k*nucl_spher_pos[0]*dfdf;
   }
   else if(basis_func_type=="2px")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*(std::complex<double>(0,-1)*sin(phi)+k*nucl_spher_pos[0]*cos(phi)*dfdf)/(k*k);
   }
   else if(basis_func_type=="2py")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*(std::complex<double>(0,1)*cos(phi)+k*nucl_spher_pos[0]*sin(phi)*dfdf)/(k*k);
   }
   else if(basis_func_type=="2pz")
   {
      return 0.55132889542179209*phase_factor*(-sqrt(contraction_zeta)*k+(2*contraction_zeta+k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*cos(thet)*dfdf/(k*k);
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.034684936146269905*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*sin(thet)*sin(thet)*(2*cos(2*phi)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(2*phi)*dfdf)/pow(k,3);
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.069369872292539811*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*sin(thet)*cos(thet)*(cos(phi)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*sin(phi)*dfdf)/pow(k,3);
   }
   else if(basis_func_type=="3d0")
   {
      return 0.010012678610436957*std::complex<double>(0,1)*nucl_spher_pos[0]*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*(6*sin(2*thet)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*(1+3*cos(2*thet))*dfdf)/pow(k,3);
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.069369872292539811*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*sin(thet)*cos(thet)*(sin(phi)+std::complex<double>(0,1)*k*nucl_spher_pos[0]*cos(phi)*dfdf)/pow(k,3);
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.034684936146269905*phase_factor*exp(-k*k/(4*contraction_zeta))*(6*contraction_zeta*k+k*k*k-6*pow(contraction_zeta,1.5)*exp(k*k/(4*contraction_zeta))*sqrt(acos(-1))*erf(k/(2*sqrt(contraction_zeta))))*sin(thet)*sin(thet)*(2*sin(2*phi)-std::complex<double>(0,1)*k*nucl_spher_pos[0]*cos(2*phi)*dfdf)/pow(k,3);
   }
   else if(basis_func_type=="4f3-")
   {
      return 0.042273611654255666*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(1+2*cos(2*phi))*sin(thet)*sin(thet)*sin(thet)*(3*std::complex<double>(0,1)*cos(3*phi)+k*nucl_spher_pos[0]*sin(3*phi)*dfdf);
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.05177438906874933*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*cos(thet)*sin(thet)*sin(thet)*(2*std::complex<double>(0,1)*cos(2*phi)+k*nucl_spher_pos[0]*sin(2*phi)*dfdf);
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.0081862496960485968*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*(3+5*cos(2*thet))*(std::complex<double>(0,1)*cos(phi)+k*nucl_spher_pos[0]*sin(phi)*dfdf);
   }
   else if(basis_func_type=="4f0")
   {
      return 0.006684044887444316*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*(-1+5*cos(2*thet))*k*nucl_spher_pos[0]*cos(thet)*dfdf;
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.0081862496960485968*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*(3+5*cos(2*thet))*(std::complex<double>(0,-1)*sin(phi)+k*nucl_spher_pos[0]*cos(phi)*dfdf);
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.05177438906874933*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*cos(thet)*sin(thet)*sin(thet)*(2*std::complex<double>(0,-1)*sin(2*phi)+k*nucl_spher_pos[0]*cos(2*phi)*dfdf);
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.042273611654255666*phase_factor*(-sqrt(contraction_zeta)*k*(30*contraction_zeta+k*k)+(60*contraction_zeta*contraction_zeta+12*contraction_zeta*k*k+k*k*k*k)*gsl_sf_dawson(k/(2*sqrt(contraction_zeta))))*sin(thet)*sin(thet)*sin(thet)*(3*std::complex<double>(0,-1)*sin(3*phi)+k*nucl_spher_pos[0]*cos(3*phi)*dfdf);
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
}
