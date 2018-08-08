#include "pice_module.hpp"
pice_set::pice_set(std::string file_address)
{

//################################INITIALIZING THE PICE PARAMETERS FROM THE HDF5 FILE #####################################
   //Initializing array size variables
   this->m_n_states_neut=new int;
   this->m_n_states_cat=new int;
   this->m_n_occ=new int;
   this->m_n_closed=new int;
   this->m_grid_size=new int;
   this->m_basis_size=new int;
   this->m_num_of_nucl=new int;
   this->m_max_cont_num=new int;
   this->m_contraction_number=new int;
      this->m_n_nucl_dim=new int;
      int *n_states_neut=new int;
      int *n_states_cat=new int;
      int *n_occ=new int;
      int *n_closed=new int;
      int *n_nucl_dim=new int;
      int *grid_size=new int;
      int *basis_size=new int;
      int *num_of_nucl=new int;
      int *max_cont_num=new int;

   //Reading array size variables
   read_output(file_address,n_states_neut,n_states_cat,n_occ,n_closed,n_nucl_dim,grid_size,num_of_nucl,basis_size);

   //Saving array size variables
   *this->m_n_states_neut=*n_states_neut;
   *this->m_n_states_cat=*n_states_cat;
   *this->m_n_occ=*n_occ;
   *this->m_n_closed=*n_closed;
   *this->m_n_nucl_dim=*n_nucl_dim;
   *this->m_grid_size=*grid_size;
   *this->m_basis_size=*basis_size;
   *this->m_num_of_nucl=*num_of_nucl;

   //Initializing contraction number vector
   int *contraction_number=new int[*this->m_basis_size];
   this->m_contraction_number=new int[*this->m_basis_size];

   //Reading contraction number vector
   read_output(file_address,n_states_neut,n_states_cat,n_occ,n_closed,n_nucl_dim,grid_size,num_of_nucl,basis_size,contraction_number);
   //Saving contraction number vector anf infering max contraction number
   for(int i=0;i!=*this->m_basis_size;i++)
   {
      this->m_contraction_number[i]=contraction_number[i];
      if(contraction_number[i] > *max_cont_num)
      {
         *max_cont_num = contraction_number[i];
      }
   }
   *this->m_max_cont_num=*max_cont_num;

   //Initializing arrays depending of variables
//   double **contraction_coeff=new double*[*basis_size];
//   double **contraction_zeta=new double*[*basis_size];
//   int *nucleus_basis_func=new int[*basis_size];
//   std::string *basis_func_type=new std::string[*basis_size];
//   double *nucl_coord=new double[*grid_size];
   this->m_nucl_coord=new double[*grid_size];
   this->m_nucl_basis_func=new int[*basis_size];
   this->m_basis_func_type=new std::string[*basis_size];
//   double ***nucl_spher_pos=new double**[*grid_size];
   double ***mo_dipole=new double**[*grid_size];
//   double **MO_coeff_neutral=new double*[*grid_size];
//   double **dyson_mo_coeff=new double*[*grid_size];
   this->m_contraction_coeff=new double*[*basis_size];
   this->m_contraction_zeta=new double*[*basis_size];
   this->m_nucl_spher_pos=new double**[*grid_size];
   this->m_mo_dipole=new double**[*grid_size];
   this->m_MO_coeff_neutral=new double*[*grid_size];
   this->m_dyson_mo_coeff=new double*[*grid_size];
   for(int i=0;i!=*basis_size;i++)
   {
//      contraction_coeff[i]=new double[*max_cont_num];
//      contraction_zeta[i]=new double[*max_cont_num];
      this->m_contraction_coeff[i]=new double[*max_cont_num];
      this->m_contraction_zeta[i]=new double[*max_cont_num];
   }

   for(int i=0;i!=*this->m_grid_size;i++)
   {
//      nucl_spher_pos[i]=new double*[*num_of_nucl];
      this->m_mo_dipole[i]=new double*[3];
      this->m_mo_dipole[i][0]=new double[(*n_occ+*n_closed)*(*n_occ+*n_closed)];
      this->m_mo_dipole[i][1]=new double[(*n_occ+*n_closed)*(*n_occ+*n_closed)];
      this->m_mo_dipole[i][2]=new double[(*n_occ+*n_closed)*(*n_occ+*n_closed)];
//      MO_coeff_neutral[i]=new double[*basis_size*(*n_occ+*n_closed)];
//      dyson_mo_coeff[i]=new double[(*n_occ+*n_closed)**n_states_neut**n_states_cat];
      this->m_nucl_spher_pos[i]=new double*[*num_of_nucl];
      this->m_MO_coeff_neutral[i]=new double[*basis_size*(*n_occ+*n_closed)];
      this->m_dyson_mo_coeff[i]=new double[(*n_occ+*n_closed)**n_states_neut**n_states_cat];

      for(int j=0;j!=*num_of_nucl;j++)
      {
         this->m_nucl_spher_pos[i][j]=new double[3];
//         nucl_spher_pos[i][j]=new double[3];
      }
   }

   read_output(file_address,this->m_n_states_neut,this->m_n_states_cat,this->m_n_occ,this->m_n_closed,this->m_n_nucl_dim,this->m_grid_size,this->m_num_of_nucl,this->m_basis_size,this->m_contraction_number,this->m_nucl_coord,this->m_nucl_spher_pos,this->m_mo_dipole,this->m_MO_coeff_neutral,this->m_dyson_mo_coeff,this->m_contraction_coeff,this->m_contraction_zeta,this->m_nucl_basis_func,this->m_basis_func_type);
/*
//################################ BUILDING THE ARRAYS THAT WILL ALLOW TO COMPUTE THE PICE #####################################
      this->k_part_s=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_p=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_d=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_f=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_s_gk=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_p_gk=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_d_gk=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_f_gk=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_s_gt=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_p_gt=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_d_gt=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_f_gt=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_s_gf=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_p_gf=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_d_gf=new std::complex<double>[this->m_contraction_zeta][nk];
      this->k_part_f_gf=new std::complex<double>[this->m_contraction_zeta][nk];
      for(int i=0;i!=;i++)
      {
         k=k_modulus[i];
         this->k_part_s[i]=(exp(-k*k/(4*this->m_contraction_zeta))/(2*pow(this->m_contraction_zeta,1.5)));
         this->k_part_p[i]=std::complex<double>(0,1)*sqrt(2.)*(k-(1./3.)*(k+pow(k,3)/(2*this->contraction_zeta))*gsl_sf_hyperg_1F1(1,2.5,-k*k/(4*this->contraction_zeta)))/(4.*sqrt(acos(-1))*pow(this->contraction_zeta,2));
      }
      */
}

bool pice_set::fill_pice(std::complex<double> *pice_x,std::complex<double> *pice_y,std::complex<double> *pice_z,int grid_index,int neut_st_index,int cat_st_index,double thet,double phi,double kp,double *ppot_vec)
{

   double *spher_pot_vec=new double[3];
   double thetp(0);
   double phip(0);
   double kpp(0);
   double *pot_vec=new double[3];
   *pice_x=0;
   *pice_y=0;
   *pice_z=0;
   if(ppot_vec==NULL)
   {
      pot_vec[0]=0;
      pot_vec[1]=0;
      pot_vec[2]=0;
   }
   else
   {
      pot_vec[0]=ppot_vec[0];
      pot_vec[1]=ppot_vec[1];
      pot_vec[2]=ppot_vec[2];
   }
   spher_pot_vec[0]=sqrt(pow(pot_vec[0],2)+pow(pot_vec[1],2)+pow(pot_vec[2],2));
   if(spher_pot_vec[0]==0)
   {
      spher_pot_vec[1]=acos(-1)/2;
      spher_pot_vec[2]=0;
   }
   else
   {
      spher_pot_vec[1]=acos(pot_vec[2]/spher_pot_vec[0]);
      if( pot_vec[0] == 0 && pot_vec[1] >= 0 )
      {
         spher_pot_vec[2]=acos(-1)/2;
      }
      else if(pot_vec[0] == 0 && pot_vec[1] < 0)
      {
         spher_pot_vec[2]=3.*acos(-1)/2; 
      }
      else if(pot_vec[0] < 0 && pot_vec[1] == 0)
      {
         spher_pot_vec[2]=acos(-1); 
      }
      else
      {
         spher_pot_vec[2]=atan(pot_vec[1]/pot_vec[0]);
      }
   }

   kpp=sqrt(pow(spher_pot_vec[0],2)+pow(kp,2)-2*spher_pot_vec[0]*kp*(sin(thet)*sin(spher_pot_vec[1])*cos(phi-spher_pot_vec[2])+cos(thet)*cos(spher_pot_vec[1])));
   thetp=acos((kp*cos(thet)-pot_vec[2])/kpp);
   phip=atan((kp*sin(thet)*sin(phi)-pot_vec[1])/(kp*sin(thet)*cos(phi)-pot_vec[0]));


   double stp(sin(thetp));
   double ctp(cos(phip));
   double spp(sin(phip));
   double cpp(cos(phip));

   for( int i=0;i!=*this->m_n_occ;i++)
   {
//      std::cout<<"here, "<<this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]<<std::endl;
       *pice_x-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
         *(
                  stp*cpp*MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  +ctp*cpp*MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
               -spp*MO_Fourier_transform_grad(i,2,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));

      *pice_y-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
               *(
                  stp*spp*MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  +ctp*spp*MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
               -cpp*MO_Fourier_transform_grad(i,2,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));

      *pice_z-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
               *(
                  ctp*MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  -stp*MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_basis_func_type,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));
          
      /*          for(int j=0;j!=n_occ;j++)
                {
                
                   temp-=dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *MO_Fourier_transform(j,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)*(mo_dipole[0][i*n_occ+j]*efield[0]+mo_dipole[1][i*n_occ+j]*efield[1]+mo_dipole[2][i*n_occ+j]*efield[2]);
                 }*/
 //            std::cout<<temp<<std::endl;
 //            
   }
}
