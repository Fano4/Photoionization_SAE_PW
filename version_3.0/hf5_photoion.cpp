#include "hf5_photoion.hpp"

bool write_output(std::string h5filename,int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int n_nucl_dim,int *grid_size,double *nucl_coord,int num_of_nucl,double **nucl_spher_pos,double **mo_dipoles_mat,int* basis_size,double *MO_coeff_neutral,double *dyson_mo_coeff,int *contraction_number,double **contraction_coeff,double **contraction_zeta,int* nucleus_basis_func,std::string *basis_func_type)
{
   using namespace H5;
      double nucl_cart_coord[*grid_size][num_of_nucl][3];
      double mo_dipole_array[*grid_size][3][*n_occ+*n_closed][*n_occ+*n_closed];
      double lcao_mo_coeff[*grid_size][*n_occ+*n_closed][*basis_size];
      double dyson_mo_coeff_array[*grid_size][*n_states_neut][*n_states_cat][*n_occ+*n_closed];
      int basis_func_type_array[*basis_size][2];
      nucl_cart_coord[0][0][0]=nucl_spher_pos[0][0]*sin(nucl_spher_pos[0][1])*cos(nucl_spher_pos[0][2]);
      nucl_cart_coord[0][0][1]=nucl_spher_pos[0][0]*sin(nucl_spher_pos[0][1])*sin(nucl_spher_pos[0][2]);
      nucl_cart_coord[0][0][2]=nucl_spher_pos[0][0]*cos(nucl_spher_pos[0][1]);
      nucl_cart_coord[0][1][0]=nucl_spher_pos[1][0]*sin(nucl_spher_pos[1][1])*cos(nucl_spher_pos[1][2]);
      nucl_cart_coord[0][1][1]=nucl_spher_pos[1][0]*sin(nucl_spher_pos[1][1])*sin(nucl_spher_pos[1][2]);
      nucl_cart_coord[0][1][2]=nucl_spher_pos[1][0]*cos(nucl_spher_pos[1][1]);
      for(int j=0;j!=*basis_size;j++)
      {
         basis_func_type_array[j][0]=spherical_harmonics_translator(basis_func_type[j],0);
         basis_func_type_array[j][1]=spherical_harmonics_translator(basis_func_type[j],1);
      }
      for(int i=0;i!=*n_occ+*n_closed;i++)
      {
         for(int j=0;j!=*n_occ+*n_closed;j++)
         {
            mo_dipole_array[0][0][i][j]=mo_dipoles_mat[0][i*(*n_occ+*n_closed)+j];
            mo_dipole_array[0][1][i][j]=mo_dipoles_mat[1][i*(*n_occ+*n_closed)+j];
            mo_dipole_array[0][2][i][j]=mo_dipoles_mat[2][i*(*n_occ+*n_closed)+j];
         }
         for(int j=0;j!=*basis_size;j++)
         {
            lcao_mo_coeff[0][i][j]=MO_coeff_neutral[i**basis_size+j];
         }
         for(int j=0;j!=*n_states_neut;j++)
         {
            for(int k=0;k!=*n_states_cat;k++)
            {
               dyson_mo_coeff_array[0][j][k][i]=dyson_mo_coeff[0**n_states_neut**n_states_cat*(*n_occ+*n_closed)+j**n_states_cat*(*n_occ+*n_closed)+k*(*n_occ+*n_closed)+i];
            }
         }
      }
   try
   {
      Exception::dontPrint();
      hsize_t dim1 (1);
      hsize_t temp(0);
      H5File file(h5filename,H5F_ACC_TRUNC);
   // CREATE ELECTRONIC STRUCT PARAMETERS GROUP
      Group electronic_struct_param(file.createGroup("/electronic_struct_param"));
      DataSpace *dataspace;
      DataSet *dataset;
      //DEBOGAGE
      /*//
      delete dim;
      dim=new hsize_t [2];
      dim[0]=*n_states_neut;
      dim[1]=2;
      double that[*n_states_neut][2];

      that[0][0]=1.5;
      that[1][0]=3.0;
      that[0][1]=2.5;
      that[1][1]=35.0;
      dataspace=new DataSpace(2,dim);
      dataset = new DataSet(electronic_struct_param.createDataSet("shittydata",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(that,PredType::NATIVE_DOUBLE);

      delete dataspace;
      delete dataset;
      */
   //create n_states_neutral_sym_dataset in electronic_struc_parameters group
      dim1=1;
      dataspace=new DataSpace(1,&dim1);
      dataset = new DataSet(electronic_struct_param.createDataSet("n_states_neut",PredType::NATIVE_INT, *dataspace));
      dataset->write(n_states_neut,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create n_states_cat_sym_dataset in electronic_struc_parameters group
      dim1=1;
      dataspace=new DataSpace(1,&dim1);
      dataset = new DataSet(electronic_struct_param.createDataSet("n_states_cat",PredType::NATIVE_INT, *dataspace));
      dataset->write(n_states_cat,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create n_closed_dataset in electronic_struc_parameters group
      dim1=1;
      dataspace=new DataSpace(1,&dim1);
      dataset = new DataSet(electronic_struct_param.createDataSet("n_mo_closed",PredType::NATIVE_INT, *dataspace));
      dataset->write(n_closed,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create n_occ_dataset in electronic_struc_parameters group
      dim1=1;
      dataspace=new DataSpace(1,&dim1);
      dataset = new DataSet(electronic_struct_param.createDataSet("n_mo_occ",PredType::NATIVE_INT, *dataspace));
      dataset->write(n_occ,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   // CLOSE ELECTRONIC STRUCT PARAMETERS GROUP
      electronic_struct_param.close();

   //CREATE NUCLEAR COORDINATES GROUP
      Group nuclear_coord(file.createGroup("/nuclear_coord"));
   //create n_nucl_dim in nuclear_coord group
      dim1=1;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(nuclear_coord.createDataSet("n_nucl_dim",PredType::NATIVE_INT, *dataspace));
      dataset->write(&n_nucl_dim,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create grid_size in nuclear_coord group
      dim1=1;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(nuclear_coord.createDataSet("grid_size",PredType::NATIVE_INT, *dataspace));
      dataset->write(grid_size,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create nucl_coord_array in nuclear_coord group
      temp=*grid_size;
      dim1=temp;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(nuclear_coord.createDataSet("nucl_coord_array",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(nucl_coord,PredType::NATIVE_DOUBLE);
      
      delete dataspace;
      delete dataset;
   //create nucl_coord_array in nuclear_coord group
      dim1=1;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(nuclear_coord.createDataSet("num_of_nucl",PredType::NATIVE_INT, *dataspace));
      dataset->write(&num_of_nucl,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;
   //create nucl_cartesian_coordinates in nuclear_coord group
      hsize_t dim3[3];
      temp=num_of_nucl;
      dim3[0]=temp;
      dim3[1]=3;
      temp=*grid_size;
      dim3[2]=temp;
      dataspace = new DataSpace(3,dim3);
      dataset = new DataSet(nuclear_coord.createDataSet("nucl_cartesian_array",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(nucl_cart_coord,PredType::NATIVE_DOUBLE);
      //DEBOGAGE
      /*
      for(int i=0;i!=num_of_nucl;i++)
      {
         for(int j=0;j!=3;j++)
         {
            std::cout<<"=> "<<nucl_cart_coord[0][i][j]<<std::endl;
         }
      }*/
      
      delete dataspace;
      delete dataset;
   //CLOSE NUCLEAR COORDINATES GROUP
      nuclear_coord.close();
      
   //CREATE MO DIPOLES GROUP
      Group mo_dipoles(file.createGroup("/mo_dipoles"));
   //create n_nucl_dim in mo_dipoles group
      hsize_t dim4[4];
      dim4[0]=*grid_size;
      dim4[1]=3;
      dim4[2]=(*n_closed+*n_occ);
      dim4[3]=(*n_closed+*n_occ);

      dataspace = new DataSpace(4,dim4);
      dataset = new DataSet(mo_dipoles.createDataSet("mo_dipoles_mat",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(mo_dipole_array,PredType::NATIVE_DOUBLE);
      
      delete dataspace;
      delete dataset;

   //CLOSE MO DIPOLES GROUP
      mo_dipoles.close();

   //CREATE LCAO COEFF GROUP
      Group lcao_coeff(file.createGroup("/lcao_coeff"));

      //create lcao mo coeffients of the neutral
      dim3[0]=*grid_size;
      dim3[1]=(*n_closed+*n_occ);
      dim3[2]=*basis_size;

      dataspace = new DataSpace(3,dim3);
      dataset = new DataSet(lcao_coeff.createDataSet("lcao_mo_coeff",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(lcao_mo_coeff,PredType::NATIVE_DOUBLE);
      
      delete dataspace;
      delete dataset;

      //create Dyson mo coeffients of the neutral
      dim4[0]=*grid_size;
      dim4[1]=*n_states_neut;
      dim4[2]=*n_states_cat;
      dim4[3]=*n_occ;

      dataspace = new DataSpace(4,dim4);
      dataset = new DataSet(lcao_coeff.createDataSet("dyson_mo_coeff",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(dyson_mo_coeff,PredType::NATIVE_DOUBLE);
      
      delete dataspace;
      delete dataset;
   //CLOSE LCAO COEFF GROUP
      lcao_coeff.close();
   //CREATE BASIS SET INFO GROUP
      Group basis_set_info(file.createGroup("/basis_set_info"));
      //create basis size dataset
      dim1=1;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(basis_set_info.createDataSet("basis_size",PredType::NATIVE_INT, *dataspace));
      dataset->write(basis_size,PredType::NATIVE_INT);
      
      delete dataspace;
      delete dataset;

      //create contraction number dataset
      dim1=*basis_size;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(basis_set_info.createDataSet("contraction_number",PredType::NATIVE_INT, *dataspace));
      dataset->write(contraction_number,PredType::NATIVE_INT); 
      delete dataspace;
      delete dataset;


      //Get the maximum number of contractions in the basis set to build the tables of contraction zeta and contractions coefficients
      int max_cont_num(0);
      for(int i=0;i!=*basis_size;i++)
      {
//         std::cout<<i<<"/"<<*basis_size<<","<<contraction_number[i]<<std::endl;
         if(contraction_number[i] > max_cont_num)
            max_cont_num=contraction_number[i];
      }
      //create contraction_coeff_dataset
      hsize_t dim2[2];
      dim2[0]=*basis_size;
      dim2[1]=max_cont_num;
      dataspace = new DataSpace(2,dim2);
      dataset = new DataSet(basis_set_info.createDataSet("contraction_coeff",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(contraction_coeff,PredType::NATIVE_DOUBLE); 
      delete dataspace;
      delete dataset;
      //create contraction_zeta_dataset
      dim2[0]=*basis_size;
      dim2[1]=max_cont_num;
      dataspace = new DataSpace(2,dim2);
      dataset = new DataSet(basis_set_info.createDataSet("contraction_zeta",PredType::NATIVE_DOUBLE, *dataspace));
      dataset->write(contraction_zeta,PredType::NATIVE_DOUBLE); 
      delete dataspace;
      delete dataset;
      //create nucleus_basis_function
      dim1=*basis_size;
      dataspace = new DataSpace(1,&dim1);
      dataset = new DataSet(basis_set_info.createDataSet("nucleus_basis_func",PredType::NATIVE_INT, *dataspace));
      dataset->write(nucleus_basis_func,PredType::NATIVE_INT); 
      delete dataspace;
      delete dataset;
      //create basis_function_type
      dim2[0]=*basis_size;
      dim2[1]=2;
      dataspace = new DataSpace(2,dim2);
      dataset = new DataSet(basis_set_info.createDataSet("basis_func_type",PredType::NATIVE_INT, *dataspace));
      dataset->write(basis_func_type_array,PredType::NATIVE_INT); 
      delete dataspace;
      delete dataset;
   //CLOSE BASIS SET INFO GROUP
      basis_set_info.close();

      file.close();
   }  
   catch(FileIException error)
    {
       error.printError();
       return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
       error.printError();
       return -1;
    }

    // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error)
    {
       error.printError();
       return -1;
    }
}
bool read_output(std::string h5filename,int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int *n_nucl_dim,int *grid_size,double *nucl_coord,int *num_of_nucl,double **nucl_spher_pos,double **mo_dipoles_mat,int* basis_size,double *MO_coeff_neutral,double *dyson_mo_coeff,int *contraction_number,double **contraction_coeff,double **contraction_zeta,int* nucleus_basis_func,std::string *basis_func_type)
{
   using namespace H5;
   int tester(0);
   try
   {
      Exception::dontPrint();

      H5File file(h5filename,H5F_ACC_RDONLY);

      DataSet *dataset=new DataSet(file.openDataSet("/electronic_struct_param/n_states_neut"));
      dataset->read(n_states_neut,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/electronic_struct_param/n_states_cat"));
      dataset->read(n_states_cat,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/electronic_struct_param/n_mo_occ"));
      dataset->read(n_occ,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/electronic_struct_param/n_mo_closed"));
      dataset->read(n_closed,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/nuclear_coord/n_nucl_dim"));
      dataset->read(n_nucl_dim,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/nuclear_coord/grid_size"));
      dataset->read(grid_size,PredType::NATIVE_INT);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/nuclear_coord/nucl_coord_array"));
      dataset->read(nucl_coord,PredType::NATIVE_DOUBLE);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/nuclear_coord/num_of_nucl"));
      dataset->read(num_of_nucl,PredType::NATIVE_INT);
      delete dataset;

      double nucl_cart_coord[*grid_size][*num_of_nucl][3];
      dataset=new DataSet(file.openDataSet("/nuclear_coord/nucl_cartesian_array"));
      dataset->read(nucl_cart_coord,PredType::NATIVE_DOUBLE);
      delete dataset;


      double mo_dipole_array[*grid_size][3][*n_occ+*n_closed][*n_occ+*n_closed];
      dataset=new DataSet(file.openDataSet("/mo_dipoles/mo_dipoles_mat"));
      dataset->read(mo_dipole_array,PredType::NATIVE_DOUBLE);
      delete dataset;

      double lcao_mo_coeff[*grid_size][*n_occ+*n_closed][*basis_size];
      dataset=new DataSet(file.openDataSet("/lcao_coeff/lcao_mo_coeff"));
      dataset->read(lcao_mo_coeff,PredType::NATIVE_DOUBLE);
      delete dataset;

      double dyson_mo_coeff_array[*grid_size][*n_states_neut][*n_states_cat][*n_occ+*n_closed];
      dataset=new DataSet(file.openDataSet("/lcao_coeff/dyson_mo_coeff"));
      dataset->read(dyson_mo_coeff_array,PredType::NATIVE_DOUBLE);
      delete dataset;

 //     std::cout<<*basis_size<<"!!!"<<std::endl;
      dataset=new DataSet(file.openDataSet("/basis_set_info/basis_size"));
      dataset->read(basis_size,PredType::NATIVE_INT);
      delete dataset;
 //     std::cout<<*basis_size<<"!!!"<<std::endl;

      dataset=new DataSet(file.openDataSet("/basis_set_info/contraction_number"));
      dataset->read(contraction_number,PredType::NATIVE_INT);
      delete dataset;

      int max_cont_num(0);
      for(int i=0;i!=*basis_size;i++)
      {
         if(contraction_number[i] > max_cont_num)
            max_cont_num=contraction_number[i];
      }

      dataset=new DataSet(file.openDataSet("/basis_set_info/contraction_coeff"));
      dataset->read(contraction_coeff,PredType::NATIVE_DOUBLE);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/basis_set_info/contraction_zeta"));
      dataset->read(contraction_zeta,PredType::NATIVE_DOUBLE);
      delete dataset;

      dataset=new DataSet(file.openDataSet("/basis_set_info/nucleus_basis_func"));
      dataset->read(nucleus_basis_func,PredType::NATIVE_INT);
      delete dataset;

      int basis_func_type_array[*basis_size][2];
      dataset=new DataSet(file.openDataSet("/basis_set_info/basis_func_type"));
      dataset->read(basis_func_type_array,PredType::NATIVE_INT);
      delete dataset;

      for(int j=0;j!=*basis_size;j++)
      {
         basis_func_type[j]=inverse_spherical_harmonics_translator(basis_func_type_array[j][0],basis_func_type_array[j][1]);
      }
      for(int i=0;i!=*n_occ+*n_closed;i++)
      {
         for(int j=0;j!=*n_occ+*n_closed;j++)
         {
            mo_dipoles_mat[0][i*(*n_occ+*n_closed)+j]=mo_dipole_array[0][0][i][j];
            mo_dipoles_mat[1][i*(*n_occ+*n_closed)+j]=mo_dipole_array[0][1][i][j];
            mo_dipoles_mat[2][i*(*n_occ+*n_closed)+j]=mo_dipole_array[0][2][i][j];
         }
         for(int j=0;j!=*basis_size;j++)
         {
            MO_coeff_neutral[i**basis_size+j]=lcao_mo_coeff[0][i][j];
         }
         for(int j=0;j!=*n_states_neut;j++)
         {
            for(int k=0;k!=*n_states_cat;k++)
            {
               dyson_mo_coeff[0**n_states_neut**n_states_cat*(*n_occ+*n_closed)+j**n_states_cat*(*n_occ+*n_closed)+k*(*n_occ+*n_closed)+i]=dyson_mo_coeff_array[0][j][k][i];
            }
         }
      }
      for(int i=0;i!=*num_of_nucl;i++)
      {
         nucl_spher_pos[i][0]=sqrt(pow(nucl_cart_coord[0][i][0],2)+pow(nucl_cart_coord[0][i][1],2)+pow(nucl_cart_coord[0][i][2],2));
       if(nucl_spher_pos[i][0]==0)
       {
          nucl_spher_pos[i][1]=0;
          nucl_spher_pos[i][2]=0;
       }
       else if(nucl_cart_coord[0][i][0]==0)
       {
          nucl_spher_pos[i][1]=acos(nucl_cart_coord[0][i][2]/nucl_spher_pos[i][0]);
          if(nucl_cart_coord[0][i][1]==0)
          {
             nucl_spher_pos[i][2]=0;
          }
          else if (nucl_cart_coord[0][i][1]>0)
          {
             nucl_spher_pos[i][2]=acos(-1)/2;
          }
          else
          {
             nucl_spher_pos[i][2]=3*acos(-1)/2;
          }
       }
       else
       {
           nucl_spher_pos[i][1]=acos(nucl_cart_coord[0][i][2]/nucl_spher_pos[i][0]);
           nucl_spher_pos[i][2]=atan(nucl_cart_coord[0][i][1]/nucl_cart_coord[0][i][0]);
       }
      }
   }
   catch(FileIException error)
    {
       error.printError();
       return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
       error.printError();
       return -1;
    }

    // catch failure caused by the DataSpace operations
    catch(DataSpaceIException error)
    {
       error.printError();
       return -1;
    }
   
}
int spherical_harmonics_translator(std::string basis_func_type,bool component)
{
   if(basis_func_type=="1s")
   {
      if(!component)
         return 0;
      else
         return 0;
   }
   else if(basis_func_type=="2px")
   {
      if(!component)
         return 1;
      else
         return 1;
   }
   else if(basis_func_type=="2py")
   {
      if(!component)
         return 1;
      else
         return -1;
   }
   else if(basis_func_type=="2pz")
   {
      if(!component)
         return 1;
      else
         return 0;
   }
   else if(basis_func_type=="3d2-")
   {
      if(!component)
         return 2;
      else
         return -2;
   }
   else if(basis_func_type=="3d1-")
   {
      if(!component)
         return 2;
      else
         return -1;
   }
   else if(basis_func_type=="3d0")
   {
      if(!component)
         return 2;
      else
         return 0;
   }
   else if(basis_func_type=="3d1+")
   {
      if(!component)
         return 2;
      else
         return 1;
   }
   else if(basis_func_type=="3d2+")
   {
      if(!component)
         return 2;
      else
         return 2;
   }
   else if(basis_func_type=="4f3-")
   {
      if(!component)
         return 3;
      else
         return -3;
   }
   else if(basis_func_type=="4f2-")
   {
      if(!component)
         return 3;
      else
         return -2;
   }
   else if(basis_func_type=="4f1-")
   {
      if(!component)
         return 3;
      else
         return -1;
   }
   else if(basis_func_type=="4f0")
   {
      if(!component)
         return 3;
      else
         return 0;
   }
   else if(basis_func_type=="4f1+")
   {
      if(!component)
         return 3;
      else
         return 1;
   }
   else if(basis_func_type=="4f2+")
   {
      if(!component)
         return 3;
      else
         return 2;
   }
   else if(basis_func_type=="4f3+")
   {
      if(!component)
         return 3;
      else
         return 3;
   }
   else
   {
      std::cout<<"Spherical harmonics not recognized in basis function : "<<basis_func_type.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }

}
std::string inverse_spherical_harmonics_translator(int l,int m)
{
   switch (l)
   {
      case 0:
         return "1s";
         break;
      case 1:
         if(m==-1)
            return "2py";
         else if(m==0)
            return "2pz";
         else if(m==1)
            return "2px";
         else 
            std::cout<<"Spherical harmonics not recognized in basis function : l = "<<l<<", m = "<<m<<std::endl;
            exit(EXIT_FAILURE);
      case 2:
         if(m==-2)
            return "3d2-";
         else if(m==-1)
            return "3d1-";
         else if(m==0)
            return "3d0";
         else if(m==1)
            return "3d1+";
         else if(m==2)
            return "3d2+";
         else 
            std::cout<<"Spherical harmonics not recognized in basis function : l = "<<l<<", m = "<<m<<std::endl;
            exit(EXIT_FAILURE);
      case 3:
         if(m==-3)
            return "4f3-";
         else if(m==-2)
            return "4f2-";
         else if(m==-1)
            return "4f1-";
         else if(m==0)
            return "4f0";
         else if(m==1)
            return "4f1+";
         else if(m==2)
            return "4f2+";
         else if(m==3)
            return "4f3+";
         else
            std::cout<<"Spherical harmonics not recognized in basis function : l = "<<l<<", m = "<<m<<std::endl;
            exit(EXIT_FAILURE);
      default:
      std::cout<<"Spherical harmonics not recognized in basis function : l = "<<l<<", m = "<<m<<std::endl;
      exit(EXIT_FAILURE);

   }

}
