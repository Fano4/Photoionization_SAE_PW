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

      *max_cont_num=0;

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

   this->m_angular_mom_numbers=new int*[*this->m_basis_size];

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
      this->m_angular_mom_numbers[i]=new int[2];
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

   for(int i=0;i!=*this->m_basis_size;i++)
   {
       this->m_angular_mom_numbers[i][0]=l_number(this->m_basis_func_type[i]); 
       this->m_angular_mom_numbers[i][1]=ml_number(this->m_basis_func_type[i],this->m_angular_mom_numbers[i][0]); 
//      std::cout<<this->m_basis_func_type[i].c_str()<<" -> "<<this->m_angular_mom_numbers[i][0]<<","<<this->m_angular_mom_numbers[i][1]<<std::endl;
   }
   std::cout<<"READING OF HF5 FILE DONE"<<std::endl;

   delete n_states_neut;
   delete n_states_cat;
   delete n_occ;
   delete n_closed;
   delete n_nucl_dim;
   delete grid_size;
   delete basis_size;
   delete num_of_nucl;
   delete [] contraction_number;
   delete [] mo_dipole;
}

bool pice_set::fill_pice(std::complex<double> *pice_x,std::complex<double> *pice_y,std::complex<double> *pice_z,int grid_index,int neut_st_index,int cat_st_index,double thet,double phi,double kp,double *ppot_vec)
{

   double *spher_pot_vec=new double[3];
   double thetp(0);
   double phip(0);
   double kpp(0);
   std::complex<double>temp(0,0);

   double xp(0);
   double yp(0);
   double zp(0);


 //  std::cout<<" !!!!! ANGLES AND RADII !!!! "<<kp<<" , "<<thet<<" , "<<phi<<std::endl;

   if(ppot_vec==NULL)
   {
      xp=kp*sin(thet)*cos(phi);
      yp=kp*sin(thet)*sin(phi);
      zp=kp*cos(thet);
   }
   else
   {

      xp=kp*sin(thet)*cos(phi)-ppot_vec[0];
      yp=kp*sin(thet)*sin(phi)-ppot_vec[1];
      zp=kp*cos(thet)-ppot_vec[2];
   }
   kpp=sqrt(pow(xp,2)+pow(yp,2)+pow(zp,2));
   if(kpp==0)
   {
      thetp=0;
      phip=0;
   }
   else
   {
      thetp=acos(zp/kpp);
      if( xp == 0 && yp >= 0 )
      {
         phip=acos(-1)/2;
      }
      else if(xp == 0 && yp < 0)
      {
         phip=3.*acos(-1)/2; 
      }
      else if(xp < 0 && yp == 0)
      {
         phip=acos(-1); 
      }
      else
      {
         phip=atan2(yp,xp);
         if(phip<0)
            phip+=2*acos(-1);
      }
   }
   delete []spher_pot_vec;

   double stp(sin(thetp));
   double ctp(cos(thetp));
   double spp(sin(phip));
   double cpp(cos(phip));

   *pice_x=0;
   *pice_y=0;
   *pice_z=0;
   for( int i=0;i!=*this->m_n_occ;i++)
   {
 //     std::cout<<"Dyson: "<<i<<" = "<<this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]<<std::endl<<"nucl position "<<this->m_nucl_spher_pos[grid_index][0][0]<<" , "<<this->m_nucl_spher_pos[grid_index][0][1]<<","<<this->m_nucl_spher_pos[grid_index][0][2]<<" ; "<<this->m_nucl_spher_pos[grid_index][1][0]<<" , "<<this->m_nucl_spher_pos[grid_index][1][1]<<","<<this->m_nucl_spher_pos[grid_index][1][2]<<std::endl;

      if(this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i] != 0)
      {
       *pice_x-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
         *(
                  stp * cpp * MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  + ctp * cpp * MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  -spp * MO_Fourier_transform_grad(i,2,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));

      *pice_y-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
               *(
                    stp * spp * MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  + ctp * spp * MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  + cpp * MO_Fourier_transform_grad(i,2,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));

      *pice_z-=std::complex<double>(0,1)*this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
               *(
                    ctp * MO_Fourier_transform_grad(i,0,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size)
                  - stp * MO_Fourier_transform_grad(i,1,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size));
         
                for(int j=0;j!=*this->m_n_occ;j++)
                { 
                      temp=this->m_dyson_mo_coeff[grid_index][neut_st_index**this->m_n_states_cat**this->m_n_occ+cat_st_index**this->m_n_occ+i]
                         *MO_Fourier_transform(j,kpp,thetp,phip,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size);
                               *pice_x-=temp*this->m_mo_dipole[grid_index][0][i**this->m_n_occ+j];
                               *pice_y-=temp*this->m_mo_dipole[grid_index][1][i**this->m_n_occ+j];
                               *pice_z-=temp*this->m_mo_dipole[grid_index][2][i**this->m_n_occ+j];
                 }
 //            std::cout<<temp<<std::endl;
      }
   }
   return 0;
//   std::cout<<"PROBE "<<kpp<<","<<thetp<<","<<phip<<" : "<<*pice_z<<std::endl;
}
//###############################################################################
//
//
//###############################################################################
int l_number(std::string bas_func_type)
{
   if(bas_func_type == "1s")
   {
      return 0;
   }
   else if(bas_func_type=="2px" || bas_func_type=="2py" || bas_func_type=="2pz" )
   {
      return 1;
   }
   else if(bas_func_type=="3d2-" || bas_func_type=="3d1-" || bas_func_type=="3d0" || bas_func_type=="3d1+" || bas_func_type=="3d2+" )
   {
      return 2;
   }
   else if( bas_func_type=="4f3-" || bas_func_type=="4f2-" || bas_func_type=="4f1-" || bas_func_type=="4f0" || bas_func_type=="4f1+" || bas_func_type=="4f2+" || bas_func_type=="4f3+")
   {
      return 3;
   }
   else
   {
      std::cout<<"ERROR IN L VALUE DETERMINATION : "<<bas_func_type.c_str()<<std::endl;
      exit(EXIT_SUCCESS);
   }
}
int ml_number(std::string bas_func_type,int l)
{
   switch(l) {

   case 0:
      return 0;
   case 1:
      if(bas_func_type=="2px")
         return 1;
      else if(bas_func_type=="2pz")
         return 0;
      else if(bas_func_type=="2py")
         return -1;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   case 2:
      if(bas_func_type=="3d2-")
         return -2;
      else if(bas_func_type=="3d1-")
         return -1;
      else if(bas_func_type=="3d0")
         return 0;
      else if(bas_func_type=="3d1+")
         return 1;
      else if(bas_func_type=="3d2+")
         return 2;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   case 3:
      if(bas_func_type=="4f3-")
         return -3;
      else if(bas_func_type=="4f2-")
         return -2;
      else if(bas_func_type=="4f1-")
         return -1;
      else if(bas_func_type=="4f0")
         return 0;
      else if(bas_func_type=="4f1+")
         return 1;
      else if(bas_func_type=="4f2+")
         return 2;
      else if(bas_func_type=="4f3+")
         return 3;
      else 
      {
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);
      }
   default:
         std::cout<<"ERROR IN SPHERICAL HARMONICS READING in ml_number : "<<bas_func_type.c_str()<<std::endl;
         exit(EXIT_FAILURE);

   }
}
int pice_set::n_occ()
{
   return *this->m_n_occ;
}
double pice_set::mo_value(double x,double y,double z,int mo_index,int grid_index)
{
   return MO_value(mo_index,x,y,z,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size);
}
std::complex<double> pice_set::mo_ft_value(double k,double thet,double phi,int mo_index,int grid_index)
{
   return MO_Fourier_transform(mo_index,k,thet,phi,this->m_nucl_spher_pos[grid_index],this->m_nucl_basis_func,this->m_contraction_number,this->m_contraction_coeff,this->m_contraction_zeta,this->m_angular_mom_numbers,this->m_MO_coeff_neutral[grid_index],*this->m_basis_size);
}
