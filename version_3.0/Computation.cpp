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

std::complex<double> MO_Fourier_transform_grad( int mo_index,int comp, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size)
{
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT_grad(i,comp,k,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,basis_func_type);
   }
   return value;
}
std::complex<double> AO_FT_grad(int ao_index,int comp,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type)
{
   std::complex<double> value(0,0);
   switch(comp)
   {
      case 0:
       //  std::cout<<contraction_number[ao_index]<<std::endl;//DEBOGAGE
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FT_grad_k(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index]);
//           std::cout<<k<<","<<thet<<","<<phi<<","<<nucl_spher_pos[0]<<","<<contraction_zeta[ao_index][i]<<","<<basis_func_type[ao_index]<<std::endl;
       }
       break;
      case 1:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FT_grad_thet(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index]);
       }
       break;
      case 2:
       for(int i=0;i!=contraction_number[ao_index];i++)
       {
           value+=contraction_coeff[ao_index][i]*contraction_FT_grad_phi(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index]);
       }
       break;
      default:
       std::cout<<"invalid component of grad vector required"<<std::endl;
       exit(EXIT_FAILURE);
       break;
   }
   return value;

}
std::complex<double> MO_Fourier_transform( int mo_index, double k, double thet, double phi,double **nucl_spher_pos,int *nucl_basis_func,int* contraction_number,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type,double *MO_neut_basis_coeff,int basis_size)
{
   std::complex<double> value(0);
   for(int i=0;i!=basis_size;i++)
   {
      value+=MO_neut_basis_coeff[mo_index*basis_size+i]*AO_FT(i,k,thet,phi,contraction_number,nucl_spher_pos[nucl_basis_func[i]-1],contraction_coeff,contraction_zeta,basis_func_type);
//      std::cout<<nucl_basis_func[i]-1<<"!!"<<nucl_spher_pos[nucl_basis_func[i]-1][0]<<","<<nucl_spher_pos[nucl_basis_func[i]-1][1]<<","<<nucl_spher_pos[nucl_basis_func[i]-1][2]<<std::endl;
   }
//   exit(EXIT_SUCCESS);
   return value;
}
std::complex<double> AO_FT(int ao_index,double k, double thet, double phi,int *contraction_number,double *nucl_spher_pos,double **contraction_coeff,double **contraction_zeta,std::string *basis_func_type)
{
   std::complex<double> value(0);
   for(int i=0;i!=contraction_number[ao_index];i++)
   {
      value+=contraction_coeff[ao_index][i]*contraction_FT(k,thet,phi,nucl_spher_pos,contraction_zeta[ao_index][i],basis_func_type[ao_index]);
   }
   return value;

}
std::complex<double> contraction_FT( double k, double thet, double phi,double* nucl_spher_pos,double contraction_zeta,std::string basis_func_type)
{
   using namespace std;

   double temp(0);
   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(std::complex<double>(0,-1)*k*nucl_spher_pos[0]*f));
   complex<double> k_part_s(0);
   complex<double> k_part_p(0);
   complex<double> k_part_d(0);
   complex<double> k_part_f(0);
   
   if(basis_func_type=="1s")
   {
      k_part_s=(exp(-k*k/(4*contraction_zeta))/(2*pow(contraction_zeta,1.5)));
   }
   else if(basis_func_type=="2px" || basis_func_type=="2py" || basis_func_type=="2pz")
   {
      k_part_p=std::complex<double>(0,1)*sqrt(2.)*(k-(1./3.)*(k+pow(k,3)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*contraction_zeta)))/(4.*sqrt(acos(-1))*pow(contraction_zeta,2));
   }
   else if(basis_func_type=="3d2-" || basis_func_type=="3d1-" || basis_func_type=="3d0" || basis_func_type=="3d1+" || basis_func_type=="3d2+")
   {
      k_part_d=-sqrt(2.)*(pow(k,2)*exp(-pow(k,2)/(4.*contraction_zeta))*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta)))/(40.*pow(contraction_zeta,2.5));
   }
   else if(basis_func_type=="3f3-" || basis_func_type=="3f2-" || basis_func_type=="3f1-" || basis_func_type=="3f0" || basis_func_type=="3f1+" || basis_func_type=="3f2+" || basis_func_type=="3f3+")
   {
      k_part_f=std::complex<double>(0,-1)*sqrt(2.)*(pow(k,3)/(24.*pow(contraction_zeta,3))+k/(4.*pow(contraction_zeta,2))-(pow(k,5)/(240.*pow(contraction_zeta,4))+pow(k,3)/(20.*pow(contraction_zeta,3))+k/(4.*pow(contraction_zeta,2)))*gsl_sf_hyperg_1F1(1,3.5,-k*k/(4*contraction_zeta)))/(sqrt(acos(-1)));
   }
   


  /* 
      std::cout<<"PROOOOOOOOOOOOOOOOBE#########################"<<std::endl;
      std::complex<double> value;
      contraction_zeta=0.05;
      nucl_spher_pos[0]=0;
   //for(int i=0;i!=1024;i++)
   //{
    //  k=1e-05+i*2.71/1024;
//      thet=acos(-1);
      //k=2.71*(i+1.)/1024;
//      phi=acos(-1)/2;
   //   phase_factor=sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]);
   //   dfdt=cos(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])-sin(thet)*cos(nucl_spher_pos[1]);
      //dfdf=(-sin(thet)*sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
      //dfdf2=(-sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
             //std::cout<<k*k*27.211/2<<","<<real(value)<<","<<imag(value)<<std::endl;
      for(int i=0;i!=100;i++)
      {
         for(int j=0;j!=50;j++)
         {
            thet=acos(-1)*i/100;
            phi=2*acos(-1)*j/50;

            value=0.5*sqrt(15./acos(-1))*sin(thet)*sin(2.*phi);
             std::cout<<thet<<"    "<<phi<<"    "<<real(value)<<"    "<<imag(value)<<std::endl;
         }std::cout<<std::endl;
   }
   exit(EXIT_SUCCESS);

*/
   if(basis_func_type=="1s")
   {
      return 0.5*(1/sqrt(acos(-1)))*phase_factor*k_part_s;
   }
   else if(basis_func_type=="2px")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*sin(thet)*cos(phi)*k_part_p;
   }
   else if(basis_func_type=="2py")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*sin(thet)*sin(phi)*k_part_p;
   }
   else if(basis_func_type=="2pz")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*cos(thet)*k_part_p;
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.25*sqrt(15./acos(-1))*pow(sin(thet),2)*sin(2.*phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.25*sqrt(15./acos(-1))*sin(2*thet)*sin(phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d0")
   {
      return 0.25*sqrt(5./acos(-1))*(2*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.25*sqrt(15./acos(-1))*sin(2.*thet)*cos(phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.25*sqrt(15./acos(-1))*cos(2*phi)*pow(sin(thet),2)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="4f3-")
   {
      return 0.25*sqrt(35/(2*acos(-1)))*(3*pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*sin(phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*sin(2.*phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f0")
   {
      return 0.25*sqrt(7./acos(-1))*cos(thet)*(2*cos(thet)-3*pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.5*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*cos(2*phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.25*sqrt(35/(2*acos(-1)))*(pow(sin(thet),2)*pow(cos(phi),2)-3.*pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*cos(phi)*phase_factor*k_part_f;
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

   complex<double> f(sin(thet)*sin(nucl_spher_pos[1])*cos(phi-nucl_spher_pos[2])+cos(thet)*cos(nucl_spher_pos[1]));
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));
   complex<double> k_part_s(0);
   complex<double> k_part_p(0);
   complex<double> k_part_d(0);
   complex<double> k_part_f(0);
   double temp(0);

   if(basis_func_type=="1s")
   {
       k_part_s=(exp(-k*k/(4*contraction_zeta))*(-std::complex<double>(0,1)*nucl_spher_pos[0]*f/pow(2*contraction_zeta,1.5)-k/(8*sqrt(2)*pow(contraction_zeta,2.5))));
   }
   else if(basis_func_type=="2px" || basis_func_type=="2py" || basis_func_type=="2pz")
   {
      k_part_p=sqrt(2/acos(-1))*(1./(4*pow(contraction_zeta,2)))*(nucl_spher_pos[0]*f*(k-(1./3.)*(k+pow(k,3)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*contraction_zeta)))+std::complex<double>(0,1)*(1-(1./3.)*(1+3*pow(k,2)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*contraction_zeta))+(k/(15*contraction_zeta))*(k+pow(k,3)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(2,3.5,-k*k/(4*contraction_zeta))));
   }

   else if(basis_func_type=="3d2-" || basis_func_type=="3d1-" || basis_func_type=="3d0" || basis_func_type=="3d1+" || basis_func_type=="3d2+")
   {
      k_part_d=(sqrt(2)/(40.*pow(contraction_zeta,2.5)))*(std::complex<double>(0,1)*nucl_spher_pos[0]*f*(pow(k,2)*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta)))-2*k*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta))+(2*pow(k,3)/(4*contraction_zeta))*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta))-(pow(k,3)/(7*contraction_zeta))*gsl_sf_hyperg_1F1(2,4.5,k*k/(4*contraction_zeta)))*exp(-k*k/(4*contraction_zeta));
   }


   else if(basis_func_type=="3f3-" || basis_func_type=="3f2-" || basis_func_type=="3f1-" || basis_func_type=="3f0" || basis_func_type=="3f1+" || basis_func_type=="3f2+" || basis_func_type=="3f3+")
   {
        k_part_f=std::complex<double>(0,-1)*sqrt(2./acos(-1))*(std::complex<double>(0,-1)*nucl_spher_pos[0]*f*(pow(k,3)/(24*pow(contraction_zeta,3))+k/(4*pow(contraction_zeta,2))-(pow(k,5)/(240*pow(contraction_zeta,4))+pow(k,3)/(20*pow(contraction_zeta,3))+k/(4*contraction_zeta*contraction_zeta))*gsl_sf_hyperg_1F1(1,3.5,-k*k/(4*contraction_zeta)))+(pow(k,2)/(8*pow(contraction_zeta,3))+1/(4*pow(contraction_zeta,2))-(pow(k,4)/(48*pow(contraction_zeta,4))+3*pow(k,2)/(20*pow(contraction_zeta,3))+1/(4*pow(contraction_zeta,2)))*gsl_sf_hyperg_1F1(1,3.5,-k*k/(4*contraction_zeta))+(pow(k,6)/(1680*pow(contraction_zeta,5))+pow(k,4)/(140*pow(contraction_zeta,4))+pow(k,2)/(28*pow(contraction_zeta,3)))*gsl_sf_hyperg_1F1(2,4.5,-k*k/(4*contraction_zeta))));
   }

   if(basis_func_type=="1s")
   {
      return 0.5*(1/sqrt(acos(-1)))*phase_factor*k_part_s;
   }
   else if(basis_func_type=="2px")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*sin(thet)*cos(phi)*k_part_p;
   }
   else if(basis_func_type=="2py")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*sin(thet)*sin(phi)*k_part_p;
   }
   else if(basis_func_type=="2pz")
   {
      return sqrt(3./(4.*acos(-1)))*phase_factor*cos(thet)*k_part_p;
   }
   else if(basis_func_type=="3d2-")
   {
      return 0.25*sqrt(15./acos(-1))*pow(sin(thet),2)*sin(2.*phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1-")
   {
      return 0.25*sqrt(15./acos(-1))*sin(2*thet)*sin(phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d0")
   {
      return 0.25*sqrt(5./acos(-1))*(2*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1+")
   {
      return 0.25*sqrt(15./acos(-1))*sin(2.*thet)*cos(phi)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d2+")
   {
      return 0.25*sqrt(15./acos(-1))*cos(2*phi)*pow(sin(thet),2)*phase_factor*k_part_d;
   }
   else if(basis_func_type=="4f3-")
   {
      return 0.25*sqrt(35/(2*acos(-1)))*(3*pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*sin(phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2-")
   {
      return 0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*sin(2.*phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1-")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f0")
   {
      return 0.25*sqrt(7./acos(-1))*cos(thet)*(2*pow(cos(thet),2)-3*pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1+")
   {
      return 0.25*sqrt(21./(2*acos(-1)))*sin(thet)*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2+")
   {
      return 0.5*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*cos(2*phi)*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f3+")
   {
      return 0.25*sqrt(35/(2*acos(-1)))*(pow(sin(thet),2)*pow(cos(phi),2)-3.*pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*cos(phi)*phase_factor*k_part_f;
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
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));
   complex<double> k_part_s(0);
   complex<double> k_part_p(0);
   complex<double> k_part_d(0);
   complex<double> k_part_f(0);
   double temp(0);
   if(basis_func_type=="1s")
   {
      k_part_s=(exp(-k*k/(4*contraction_zeta))/(2*pow(contraction_zeta,1.5)));
   }
   else if(basis_func_type=="2px" || basis_func_type=="2py" || basis_func_type=="2pz")
   {
      k_part_p=std::complex<double>(0,1)*sqrt(2.)*(1-(1./3.)*(1+pow(k,2)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*contraction_zeta)))/(4.*sqrt(acos(-1))*pow(contraction_zeta,2));
   }
   else if(basis_func_type=="3d2-" || basis_func_type=="3d1-" || basis_func_type=="3d0" || basis_func_type=="3d1+" || basis_func_type=="3d2+")
   {
      k_part_d=-sqrt(2.)*(k*exp(-pow(k,2)/(4.*contraction_zeta))*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta)))/(40.*pow(contraction_zeta,2.5));
   }
   else if(basis_func_type=="3f3-" || basis_func_type=="3f2-" || basis_func_type=="3f1-" || basis_func_type=="3f0" || basis_func_type=="3f1+" || basis_func_type=="3f2+" || basis_func_type=="3f3+")
   {
      k_part_f=std::complex<double>(0,-1)*sqrt(2.)*(pow(k,2)/(24.*pow(contraction_zeta,3))+1/(4.*pow(contraction_zeta,2))-(pow(k,4)/(240.*pow(contraction_zeta,4))+pow(k,2)/(20.*pow(contraction_zeta,3))+(1/4.*pow(contraction_zeta,2)))*gsl_sf_hyperg_1F1(1,3.5,-k*k/(4*contraction_zeta)))/(sqrt(acos(-1)));
   }
   if(basis_func_type=="1s")
   {
      return std::complex<double>(0,-1)*nucl_spher_pos[0]*dfdt*0.5*(1/sqrt(acos(-1)))*phase_factor*k_part_s;
   }
   else if(basis_func_type=="2px")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(sqrt(3./(4.*acos(-1)))*sin(thet)*cos(phi))+sqrt(3./(4.*acos(-1)))*cos(thet)*cos(phi))*phase_factor*k_part_p;
   }
   else if(basis_func_type=="2py")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(sqrt(3./(4.*acos(-1)))*sin(thet)*sin(phi))+sqrt(3./(4.*acos(-1)))*cos(thet)*sin(phi))*phase_factor*k_part_p;
   }
   else if(basis_func_type=="2pz")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(sqrt(3./(4.*acos(-1)))*cos(thet))-sqrt(3./(4.*acos(-1)))*sin(thet))*phase_factor*k_part_p;
   }
   else if(basis_func_type=="3d2-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(15./acos(-1))*pow(sin(thet),2)*sin(2.*phi))+0.25*sqrt(15./acos(-1))*sin(2*thet)*sin(2.*phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(15./acos(-1))*sin(2*thet)*sin(phi))+0.5*sqrt(15./acos(-1))*cos(2*thet)*sin(phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d0")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(5./acos(-1))*(2*pow(cos(thet),2)-pow(sin(thet),2)))-0.75*sqrt(5./acos(-1))*(sin(2*thet)))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(15./acos(-1))*sin(2.*thet)*cos(phi))+0.5*sqrt(15./acos(-1))*cos(2*thet)*cos(phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d2+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(15./acos(-1))*cos(2*phi)*pow(sin(thet),2))+0.25*sqrt(15./acos(-1))*sin(2*thet)*cos(2.*phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="4f3-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(35/(2*acos(-1)))*(3*pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2)))+0.75*sqrt(35/(2*acos(-1)))*pow(sin(thet),2)*cos(thet)*(3*pow(cos(phi),2)-pow(sin(phi),2))*sin(phi))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*sin(2.*phi))+0.25*sqrt(105./acos(-1))*(2.*sin(thet)*pow(cos(thet),2)-pow(sin(thet),3))*sin(2.*phi))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(21./(2*acos(-1)))*sin(thet)*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))+0.25*sqrt(21./(2*acos(-1)))*sin(phi)*(4*pow(cos(thet),3)-11.*pow(sin(thet),2)*cos(thet)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f0")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(7./acos(-1))*cos(thet)*(2*pow(cos(thet),2)-3*pow(sin(thet),2)))-0.75*sqrt(7./acos(-1))*(2*sin(2*thet)*cos(thet)-pow(sin(thet),3)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(21./(2*acos(-1)))*sin(thet)*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))+0.25*sqrt(21./(2*acos(-1)))*cos(phi)*(4*pow(cos(thet),3)-11*pow(sin(thet),2)*cos(thet)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.5*sqrt(105./acos(-1))*pow(sin(thet),2)*cos(thet)*cos(2*phi))+0.5*sqrt(105./acos(-1))*(2*sin(thet)*pow(cos(thet),2)-pow(sin(thet),3))*cos(2.*phi))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f3+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdt*(0.25*sqrt(35/(2*acos(-1)))*(pow(sin(thet),2)*pow(cos(phi),2)-3.*pow(sin(thet),2)*pow(sin(phi),2))*sin(thet)*cos(phi))+0.75*sqrt(35/(2*acos(-1)))*pow(sin(thet),2)*cos(thet)*(pow(cos(phi),2)-3*pow(sin(phi),2))*cos(phi))*phase_factor*k_part_f;
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
   complex<double> dfdf(-sin(thet)*sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
   complex<double> dfdf2(-sin(nucl_spher_pos[1])*sin(phi-nucl_spher_pos[2]));
   complex<double> phase_factor(exp(-std::complex<double>(0,1)*k*nucl_spher_pos[0]*f));
   complex<double> k_part_s(0);
   complex<double> k_part_p(0);
   complex<double> k_part_d(0);
   complex<double> k_part_f(0);
   double temp(0);
   if(basis_func_type=="1s")
   {
      k_part_s=(exp(-k*k/(4*contraction_zeta))/(2*pow(contraction_zeta,1.5)));
   }
   else if(basis_func_type=="2px" || basis_func_type=="2py" || basis_func_type=="2pz")
   {
      k_part_p=std::complex<double>(0,1)*sqrt(2.)*(1-(1./3.)*(1+pow(k,2)/(2*contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*contraction_zeta)))/(4.*sqrt(acos(-1))*pow(contraction_zeta,2));
   }
   else if(basis_func_type=="3d2-" || basis_func_type=="3d1-" || basis_func_type=="3d0" || basis_func_type=="3d1+" || basis_func_type=="3d2+")
   {
      k_part_d=-sqrt(2.)*(k*exp(-pow(k,2)/(4.*contraction_zeta))*gsl_sf_hyperg_1F1(1,3.5,k*k/(4*contraction_zeta)))/(40.*pow(contraction_zeta,2.5));
   }
   else if(basis_func_type=="3f3-" || basis_func_type=="3f2-" || basis_func_type=="3f1-" || basis_func_type=="3f0" || basis_func_type=="3f1+" || basis_func_type=="3f2+" || basis_func_type=="3f3+")
   {
      k_part_f=std::complex<double>(0,-1)*sqrt(2.)*(pow(k,2)/(24.*pow(contraction_zeta,3))+1/(4.*pow(contraction_zeta,2))-(pow(k,4)/(240.*pow(contraction_zeta,4))+pow(k,2)/(20.*pow(contraction_zeta,3))+(1/4.*pow(contraction_zeta,2)))*gsl_sf_hyperg_1F1(1,3.5,-k*k/(4*contraction_zeta)))/(sqrt(acos(-1)));
   }
   if(basis_func_type=="1s")
   {
      return std::complex<double>(0,-1)*nucl_spher_pos[0]*dfdf2*0.5*(1/sqrt(acos(-1)))*phase_factor*k_part_s;
   }
   else if(basis_func_type=="2px")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(sqrt(3./(4.*acos(-1)))*cos(phi))-sqrt(3./(4.*acos(-1)))*sin(phi))*phase_factor*k_part_p;
   }
   else if(basis_func_type=="2py")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(sqrt(3./(4.*acos(-1)))*sin(phi))+sqrt(3./(4.*acos(-1)))*cos(phi))*phase_factor*k_part_p;
   }
   else if(basis_func_type=="2pz")
   {
      return std::complex<double>(0,-1)*nucl_spher_pos[0]*dfdf2*sqrt(3./(4.*acos(-1)))*cos(thet)*phase_factor*k_part_p;
   }
   else if(basis_func_type=="3d2-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(15./acos(-1))*sin(thet)*sin(2.*phi))+0.5*sqrt(15./acos(-1))*sin(thet)*cos(2.*phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d1-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.5*sqrt(15./acos(-1))*cos(thet)*sin(phi))+0.5*sqrt(15./acos(-1))*cos(thet)*cos(phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d0")
   {
      return std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf2*(0.25*sqrt(5./acos(-1))*(2*pow(cos(thet),2)-pow(sin(thet),2))*phase_factor*k_part_d);
   }
   else if(basis_func_type=="3d1+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(15./acos(-1))*sin(2.*thet)*cos(phi))-0.5*sqrt(15./acos(-1))*cos(thet)*sin(phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="3d2+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(15./acos(-1))*cos(2*phi)*sin(thet))-0.5*sqrt(15./acos(-1))*sin(thet)*sin(2.*phi))*phase_factor*k_part_d;
   }
   else if(basis_func_type=="4f3-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(35/(2*acos(-1)))*(3*pow(sin(thet),2)*pow(cos(phi),2)-pow(sin(thet),2)*pow(sin(phi),2))*sin(phi))+0.75*sqrt(35/(2*acos(-1)))*pow(sin(thet),2)*(pow(cos(phi),3)-3*pow(sin(phi),2)*cos(phi)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(105./acos(-1))*sin(thet)*cos(thet)*sin(2.*phi))+0.5*sqrt(105./acos(-1))*sin(thet)*cos(thet)*cos(2.*phi))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1-")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(21./(2*acos(-1)))*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))+0.25*sqrt(21./(2*acos(-1)))*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f0")
   {
      return std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf2*(0.25*sqrt(7./acos(-1))*cos(thet)*(2*pow(cos(thet),2)-3*pow(sin(thet),2)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f1+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(21./(2*acos(-1)))*cos(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))-0.25*sqrt(21./(2*acos(-1)))*sin(phi)*(4*pow(cos(thet),2)-pow(sin(thet),2)))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f2+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.5*sqrt(105./acos(-1))*sin(thet)*cos(thet)*cos(2*phi))-sqrt(105./acos(-1))*sin(thet)*cos(thet)*sin(2.*phi))*phase_factor*k_part_f;
   }
   else if(basis_func_type=="4f3+")
   {
      return (std::complex<double>(0,-1)*k*nucl_spher_pos[0]*dfdf*(0.25*sqrt(35/(2*acos(-1)))*(pow(sin(thet),2)*pow(cos(phi),2)-3.*pow(sin(thet),2)*pow(sin(phi),2))*cos(phi))+0.75*sqrt(35/(2*acos(-1)))*(pow(sin(phi),3)-3*sin(phi)*pow(cos(phi),2))*pow(sin(thet),2))*phase_factor*k_part_f;
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
}
