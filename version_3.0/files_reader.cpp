
int n_states_reader(int *n_states_neut,int *n_states_cat,int *n_elec_neut,std::string file_address)
{
    std::string tmp_str;
    int position(0);
    double floating;
    std::ifstream molpro_output;
    
    if(!search(&position, file_address,0, "{casscf"))
    {
        std::cout<<"NUMBER OF ELECTRONIC STATES NOT FOUND IN MOLPRO OUTPUT FILE"<<std::endl;
        return 1;
    }
    
    molpro_output.open(file_address.c_str());
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    *n_states_neut=tmp_str.at(tmp_str.length()-2)-'0';
    *n_elec_neut=tmp_str.at(3)-'0';

    std::cout<<" neutral states " <<*n_states_neut<<std::endl;
    molpro_output.close();

    if(!search(&position, file_address,position, "{casscf"))
    {
        std::cout<<"NUMBER OF ELECTRONIC STATES NOT FOUND IN MOLPRO OUTPUT FILE"<<std::endl;
        return 1;
    }
    molpro_output.open(file_address.c_str());
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    molpro_output>>tmp_str;
    *n_states_cat=tmp_str.at(tmp_str.length()-2)-'0';
    
    std::cout<<" cat states " <<*n_states_cat<<std::endl;
    molpro_output.close();
    return 0;
}

int size_query(int* n_occ,int *n_closed,int* basis_size,std::string molpro_out_path,int n_sym=1)
{
    bool test(0);
    
    using namespace std;
    
    string tmp_str;
    int integer;
    double floating;
    int position(0);
    
    
    
        ifstream molpro_file;
        molpro_file.open(molpro_out_path.c_str());
        
        if (!molpro_file.is_open())
        {
            cout<<"ERROR DURING MOLPRO OUTPUT FILE ACQUISITION"<<endl<<molpro_out_path.c_str()<<endl;
            exit(EXIT_FAILURE);
        }
        
       if(n_sym==1)
       {
        
        molpro_file>>tmp_str;
        do
        {
            while (tmp_str!="NUMBER")
            {
                molpro_file>>tmp_str;
            }
            molpro_file>>tmp_str;
            molpro_file>>tmp_str;
            
            if (tmp_str=="CONTRACTIONS:")
            {
                molpro_file>>integer;
                *basis_size=integer;
                
                cout<<"Basis set size found!! Size= "<<*basis_size<<endl;//DEBOGAGE
                test=1;
            }
        }while(test!=1);//Search for the basis set size in molpro output
        molpro_file.close();

        if(!search(&position, molpro_out_path,position, "closed,"))
        {
            std::cout<<"NUMBER OF CLOSED ORBITALS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
            return 1;
        }
        molpro_file.open(molpro_out_path.c_str());

        molpro_file.seekg(position);
        molpro_file>>integer;
        *n_closed=integer;
        molpro_file.close();

        position=0;

        if(!search(&position, molpro_out_path,position, "closed,", 1, "occ,"))
        {
            std::cout<<"NUMBER OF OCCUPIED ORBITALS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
            return 1;
        }

        molpro_file.open(molpro_out_path.c_str());

        molpro_file.seekg(position);
        molpro_file>>integer;
        *n_occ=integer;
     }

       else
       {
        
        molpro_file>>tmp_str;
        do
        {
            while (tmp_str!="NUMBER")
            {
                molpro_file>>tmp_str;
            }
            molpro_file>>tmp_str;
            molpro_file>>tmp_str;
            
            if (tmp_str=="CONTRACTIONS:")
            {
                molpro_file>>integer;
                *basis_size=integer;
                molpro_file>>tmp_str;
                cout<<"Basis set size found!!"<<endl;//DEBOGAGE
                for(int i=0;i!=n_sym;i++)
                {
                   molpro_file>>tmp_str;
                   basis_size[i]=atoi(separateThem(tmp_str).c_str());
                   std::cout<<"sym "<<i+1<<" = "<<basis_size[i]<<std::endl;
                   molpro_file>>tmp_str;
                }
                
                test=1;
            }
        }while(test!=1);//Search for the basis set size in molpro output
        molpro_file.close();

        if(!search(&position, molpro_out_path,position, "closed,"))
        {
            std::cout<<"NUMBER OF CLOSED ORBITALS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
            return 1;
        }
        molpro_file.open(molpro_out_path.c_str());

        molpro_file.seekg(position);
        for(int i=0;i!=n_sym;i++)
        {
           molpro_file>>integer;
           n_closed[i]=integer;
           molpro_file>>tmp_str;
           std::cout<<"closed "<<n_closed[i]<<std::endl;
        }
           molpro_file.close();

        if(!search(&position, molpro_out_path,position, "occ,"))
        {
            std::cout<<"NUMBER OF OCCUPIED ORBITALS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
            return 1;
        }

        molpro_file.open(molpro_out_path.c_str());

        molpro_file.seekg(position);
        for(int i=0;i!=n_sym;i++)
        {
            molpro_file>>integer;
            n_occ[i]=integer;
            molpro_file>>tmp_str;
           std::cout<<"occ "<<n_occ[i]<<std::endl;
        }
     }

    molpro_file.close();

    return 0;
}


void dipole_MO(double **matrix,int* n_occ,int* basis_size,int* basis_size_sym,std::string molpro_out_path,int n_sym)
{
   //THE DMOLECULAR ORBITALS DIPOLE MOMENT MATRIX IS READ FROM MOLPRO OUTPUT FILE USING THE AO DIPOLE MOMENT MATRIX.
   //THIS MATRIX HAS A MORE COMPLICATED FORM THAN THE OVERLAP MATRIX BECAUSE IT IS NOT BLOCK-DIAGONAL!
   //EACH COMPONENT BELONNGS TO A DIFFERENT SYMMETRY AND THE STRUCTURE OF THE MATRIX REFLECTS THE SYMMETRY OF EACH COMPONENT OF THE VECTOR.
   
   //THE FIRST STEP IS TO READ THE AO TRANSITION DIPOLE MATRIX FROM MOLPRO OUTPUT FILE AND PUT IT IN GOOD FORM IN AN ARRAY, FOR EACH COMPONENT. IN THIS VERSION OF THE CODE, WE ASSUME THE FOLLOWING SYMMETRY RULES:
   // Z = 1 (A1); X = 2 (B1); Y = 3 (B2)
   //
   //     A1   B1   B2   A2
   //     ----------------
   //A1|  A1   B1   B2   A2
   //  |
   //B1|  B1   A1   A2   B2
   //  |  
   //B2|  B2   A2   A1   B1
   //  |
   //A2|  A2   B2   B1   A1      
   //  |
   
   using namespace std;
   int position(0);
   int total1(0);
   int total2(0);
   int total3(0);
   int total4(0);
   double *temp;
   double *temp2;
   double floating;
   string tmp_str;
   double *aodipole_x = new double [*basis_size**basis_size];
   double *aodipole_y = new double [*basis_size**basis_size];
   double *aodipole_z = new double [*basis_size**basis_size];
   int n_occ_tot(0);
    for(int k=0;k!=n_sym;k++)
    {
       n_occ_tot+=n_occ[k];
    } 
   double *MO_coeff_neutral=new double[*basis_size*n_occ_tot];
   for(int i=0;i!=*basis_size**basis_size;i++)
   {
      aodipole_x[i]=0;
      aodipole_y[i]=0;
      aodipole_z[i]=0;
   }

   ifstream input;

   //X COMPONENT RESEARCH
   // THE MATRIX STRUCTURE INCLUDES ALL B1 DIRECT PRODUCTS IN THE ORDER REPRESENTED IN THE PRODUCT TABLE.

   if(!search(&position, molpro_out_path,0, "MATRIX", 0, "DMX_MO"))
   {
      std::cout<<"ERROR: CANNOT FIND X DIPOLE MATRIX IN MOLPRO OUTPUT!"<<std::endl<<molpro_out_path.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }

   input.open(molpro_out_path.c_str());
   input.seekg(position);


   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //A1*B1
   total1=0;
   total2=basis_size_sym[0];

   for(int i=0;i!=basis_size_sym[0];i++)
   {
      for(int j=0;j!=basis_size_sym[1];j++)
      {
         input>>aodipole_x[(total1+i)**basis_size+j+total2];
      }
   }

   //B1*A1
   total1=basis_size_sym[0];
   total2=0;

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[1];i++)
   {
      for(int j=0;j!=basis_size_sym[0];j++)
      {
         input>>aodipole_x[(total1+i)**basis_size+j+total2];
      }
   }
   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //B2*A2
   total1=basis_size_sym[0]+basis_size_sym[1];
   total2=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];

   for(int i=0;i!=basis_size_sym[2];i++)
   {
      for(int j=0;j!=basis_size_sym[3];j++)
      {
         input>>aodipole_x[(total1+i)**basis_size+j+total2];
      }
   }

   //A2*B2
   total1=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];
   total2=basis_size_sym[0]+basis_size_sym[1];

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[3];i++)
   {
      for(int j=0;j!=basis_size_sym[2];j++)
      {
         input>>aodipole_x[(total1+i)**basis_size+j+total2];
      }
   }
   input.close();
   //Y COMPONENT RESEARCH
   position=0;
   if(!search(&position, molpro_out_path,position, "MATRIX", 0, "DMY_MO"))
   {
      std::cout<<"ERROR: CANNOT FIND Y DIPOLE MATRIX IN MOLPRO OUTPUT!"<<std::endl<<molpro_out_path.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }

   input.open(molpro_out_path.c_str());
   input.seekg(position);


   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //A1*B2
   total1=0;
   total2=basis_size_sym[0]+basis_size_sym[1];

   for(int i=0;i!=basis_size_sym[0];i++)
   {
      for(int j=0;j!=basis_size_sym[2];j++)
      {
         input>>aodipole_y[(total1+i)**basis_size+j+total2];
      }
   }

   //B1*A2
   total1=basis_size_sym[0];
   total2=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[1];i++)
   {
      for(int j=0;j!=basis_size_sym[3];j++)
      {
         input>>aodipole_y[(total1+i)**basis_size+j+total2];
      }
   }
   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //B2*A1

   total1=basis_size_sym[0]+basis_size_sym[1];
   total2=0;

   for(int i=0;i!=basis_size_sym[2];i++)
   {
      for(int j=0;j!=basis_size_sym[0];j++)
      {
         input>>aodipole_y[(total1+i)**basis_size+j+total2];
      }
   }

   //A2*B1
   total1=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];
   total2=basis_size_sym[0];

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[3];i++)
   {
      for(int j=0;j!=basis_size_sym[1];j++)
      {
         input>>aodipole_y[(total1+i)**basis_size+j+total2];
      }
   }
   input.close();
   //Z COMPONENT RESEARCH
   position=0;
   if(!search(&position, molpro_out_path,position, "MATRIX", 0, "DMZ_MO"))
   {
      std::cout<<"ERROR: CANNOT FIND Z DIPOLE MATRIX IN MOLPRO OUTPUT!"<<std::endl<<molpro_out_path.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }

   input.open(molpro_out_path.c_str());
   input.seekg(position);


   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //A1*A1
   total1=0;
   total2=0;

   for(int i=0;i!=basis_size_sym[0];i++)
   {
      for(int j=0;j!=basis_size_sym[0];j++)
      {
         input>>aodipole_z[(total1+i)**basis_size+j+total2];
      }
   }

   //B1*B1
   total1=basis_size_sym[0];
   total2=basis_size_sym[0];

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[1];i++)
   {
      for(int j=0;j!=basis_size_sym[1];j++)
      {
         input>>aodipole_z[(total1+i)**basis_size+j+total2];
      }
   }
   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;

   //B2*B2

   total1=basis_size_sym[0]+basis_size_sym[1];
   total2=basis_size_sym[0]+basis_size_sym[1];

   for(int i=0;i!=basis_size_sym[2];i++)
   {
      for(int j=0;j!=basis_size_sym[2];j++)
      {
         input>>aodipole_z[(total1+i)**basis_size+j+total2];
      }
   }

   //A2*A2
   total1=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];
   total2=basis_size_sym[0]+basis_size_sym[1]+basis_size_sym[2];

   input>>tmp_str;
   input>>tmp_str;
   input>>tmp_str;
   for(int i=0;i!=basis_size_sym[3];i++)
   {
      for(int j=0;j!=basis_size_sym[3];j++)
      {
         input>>aodipole_z[(total1+i)**basis_size+j+total2];
      }
   }
   input.close();

/*    for(int i=0;i!=*basis_size;i++)
    {
       for(int j=0;j!=*basis_size;j++)
       {
          std::cout<<aodipole_z[i**basis_size+j]<<"   ";
       }std::cout<<std::endl<<std::endl;

    }
    exit(EXIT_SUCCESS);
  */ //ONCE THE AO TRANSITION DIPOLE MATRICES ARE KNOWM, WE HAVE TO COMPUTE THE MOLECULAR ORBITALS TRANSITION DIPOLES USING THE LCAO COEFFICIENTS.
   //SEARCH FOR THE LCAO COEFFICIENTS
   position=0;
   if(!search(&position, molpro_out_path,position, "NATURAL", 0, "ORBITALS"))
   {
      std::cout<<"LCAO COEFFICIENTS NOT FOUND IN MOLPRO OUTPUT"<<std::endl;
      exit(EXIT_FAILURE);
   }
       input.open(molpro_out_path.c_str());
       input.seekg(position);
           for(int i=0;i!=7;i++)
           {
              input>>tmp_str;
           }
           total1=0;
           total2=0;
           for(int s=0;s!=n_sym;s++)
           {
            //cout<<"Coefficients of the neutral"<<endl;//DEBOGAGE
            for (int i=0; i!=2*basis_size_sym[s]; i++)
            {
                input>>tmp_str;
            }
            for (int j=total1; j!=n_occ[s]+total1; j++)
            {
                input>>tmp_str;
                input>>tmp_str;
                input>>tmp_str;
                for (int k=total2; k!=basis_size_sym[s]+total2; k++)
                {
                    input>>floating;
                    MO_coeff_neutral[j**basis_size+k]=floating;
                }//cout<<endl;
            }
            total2+=basis_size_sym[s];
            total1+=n_occ[s];
        }
    
        input.close();

    temp=new double[*basis_size*n_occ_tot];
    temp2=new double[*basis_size*n_occ_tot];
    for(int i=0;i!=*basis_size*n_occ_tot;i++)
    {
       temp[i]=0;
       temp2[i]=0;
    }

    for(int k=0;k!=n_sym;k++)
    {
       total1=0;
       total3=0;
       for(int i=0;i!=k;i++)
       {
          total1+=n_occ[i];
          total3+=basis_size_sym[i];
       }
       for(int kp=0;kp!=n_sym;kp++)
       {
          total2=0;
          total4=0;
          for(int i=0;i!=kp;i++)
          {
             total2+=n_occ[i];
             total4+=basis_size_sym[i];
          }
          for(int i=0;i!=n_occ[k];i++)
          {
             for(int j=0;j!=n_occ[kp];j++)
             {
                matrix[0][(i+total1)*n_occ_tot+j+total2]=aodipole_x[(i+total3)**basis_size+j+total4];
                matrix[1][(i+total1)*n_occ_tot+j+total2]=aodipole_y[(i+total3)**basis_size+j+total4];
                matrix[2][(i+total1)*n_occ_tot+j+total2]=aodipole_z[(i+total3)**basis_size+j+total4];
             }
          }
       }
    }

 /*   for(int i=0;i!=n_occ_tot;i++)
    {
       total1=0;
       for(int j=0;j!=*basis_size;j++)
       {
          std::cout<<setprecision(10)<<setw(15)<<MO_coeff_neutral[i**basis_size+j];
          total1++;
          if(total1%6==0)
             std::cout<<std::endl;
       }std::cout<<std::endl;

    }*/




/*    for(int i=0;i!=n_occ_tot;i++)
    {
       for(int j=0;j!=n_occ_tot;j++)
       {
          matrix[0][i*n_occ_tot+j]=0;
          matrix[1][i*n_occ_tot+j]=0;
          matrix[2][i*n_occ_tot+j]=0;
          for(int r=0;r!=*basis_size;r++)
          {
             for(int s=0;s!=*basis_size;s++)
             {
                matrix[0][i*n_occ_tot+j]+=MO_coeff_neutral[i**basis_size+r]*MO_coeff_neutral[j**basis_size+s]*aodipole_x[r**basis_size+s];
                matrix[1][i*n_occ_tot+j]+=MO_coeff_neutral[i**basis_size+r]*MO_coeff_neutral[j**basis_size+s]*aodipole_y[r**basis_size+s];
                matrix[2][i*n_occ_tot+j]+=MO_coeff_neutral[i**basis_size+r]*MO_coeff_neutral[j**basis_size+s]*aodipole_z[r**basis_size+s];
             }
          }
       }
    }*/
/*    std::cout<<"X dipole "<<n_occ_tot<<std::endl;
    transpose(MO_coeff_neutral, temp2, n_occ_tot, *basis_size);
    matrix_product(temp, aodipole_x, temp2, *basis_size, *basis_size, n_occ_tot);
    matrix_product(matrix[0], MO_coeff_neutral, temp, n_occ_tot, *basis_size, n_occ_tot); 
    std::cout<<"Y dipole"<<std::endl;
    matrix_product(temp, aodipole_y, temp2, *basis_size, *basis_size, n_occ_tot);
    matrix_product(matrix[1], MO_coeff_neutral, temp, n_occ_tot, *basis_size, n_occ_tot); 
    std::cout<<"Z dipole"<<std::endl;
    matrix_product(temp, aodipole_z, temp2, *basis_size, *basis_size, n_occ_tot);
    matrix_product(matrix[2], MO_coeff_neutral, temp, n_occ_tot, *basis_size, n_occ_tot);*/ 
/*    for(int i=0;i!=*basis_size;i++)
    {
       total1=0;
       for(int j=0;j!=*basis_size;j++)
       {
          std::cout<<setprecision(10)<<setw(15)<<aodipole_z[i**basis_size+j];
          total1++;
          if(total1%7==0)
             std::cout<<std::endl;
       }std::cout<<std::endl<<std::endl;

    }
    std::cout<<"#################################"<<std::endl;
    for(int i=0;i!=*basis_size;i++)
    {
       total1=0;
       for(int j=0;j!=n_occ_tot;j++)
       {
          std::cout<<setprecision(10)<<setw(15)<<temp2[i*n_occ_tot+j];
          total1++;
          if(total1%7==0)
             std::cout<<std::endl;
       }std::cout<<std::endl<<std::endl;

    }
  *//*  std::cout<<"#################################"<<std::endl;
    for(int i=0;i!=n_occ_tot;i++)
    {
       total1=0;
       for(int j=0;j!=n_occ_tot;j++)
       {
          std::cout<<setprecision(10)<<setw(15)<<matrix[2][i*n_occ_tot+j];
          total1++;
          if(total1%7==0)
             std::cout<<std::endl;
       }std::cout<<std::endl<<std::endl;

    }
    exit(EXIT_SUCCESS);
    */

}
int overlap_MO(double matrix[],int* n_occ,int* basis_size,int* basis_size_sym,std::string molpro_out_path,double* MO_coeff_neutral,int n_sym=1)
{
    bool test(0);
    bool test2(0);
    int *count=new int [n_sym];
    int total(0);
    int total2(0);
    int position(0);
    double *overlap;
    double *matrix2;
    double *MO_coeff_cation;
    double floating;
    double *temp;
    double *temp2;
    int n_occ_tot(0);

    for(int k=0;k!=n_sym;k++)
    {
       n_occ_tot+=n_occ[k];
    } 
//    std::cout<<n_occ_tot<<std::endl;
    using namespace std;
    
    string tmp_str;
    
    
    std::cout<<"basis size is "<<*basis_size<<std::endl;
    overlap=new double[*basis_size**basis_size];
    matrix2=new double [n_occ_tot*n_occ_tot];
    for(int i=0;i!=*basis_size**basis_size;i++)
    {
       overlap[i]=0;
    }
    for(int i=0;i!=n_occ_tot*n_occ_tot;i++)
    {
       matrix2[i]=0;
    }
    
    
    ifstream molpro_file;
    
    
   if(!search(&position, molpro_out_path,position, "MATRIX", 0, "S"))
   {
      std::cout<<"OVERLAP MATRIX NOT FOUND IN MOLPRO OUTPUT"<<std::endl;
      return -1;
   }
    molpro_file.open(molpro_out_path.c_str());
    if (!molpro_file.is_open())
    {
        cout<<"ERROR DURING MOLPRO OUTPUT FILE ACQUISITION"<<endl<<molpro_out_path.c_str()<<endl;
        return -1;
    }
    molpro_file.close();

    for(int k=0;k!=n_sym;k++)
    {
       molpro_file.open(molpro_out_path.c_str());
 //   std::cout<<"probe"<<std::endl;
       count[k]=0;
       molpro_file.seekg(position);
       molpro_file>>tmp_str;
       molpro_file>>tmp_str;
       molpro_file>>tmp_str;
       position=molpro_file.tellg();         

       do
       {
          molpro_file>>tmp_str;
          count[k]++;
       }while(tmp_str!="SYMMETRY" && tmp_str!="**********************************************************************************************************************************");
       if(tmp_str=="**********************************************************************************************************************************")
          test2=1;
       count[k]--;

      molpro_file.seekg(position);
  //     std::cout<<"$$$"<<std::endl;

       //basis_size_sym[k]=sqrt(count[k]);

       for (int i=total; i<basis_size_sym[k]+total; i++)
       {
          for(int j=total;j<basis_size_sym[k]+total;j++)
          {
//             std::cout<<"probei"<<std::endl;
             molpro_file>>tmp_str;
//             std::cout<<"probef"<<std::endl;
             floating=atof(tmp_str.c_str());
             overlap[i**basis_size+j]=floating;
//             std::cout<<"value "<<floating<<" assigned"<<std::endl;
//             std::cout<<basis_size_sym[k]<<";"<<position<<";"<<i<<";"<<j<<std::endl;
          }
                    // cout<<overlap[i**basis_size+total]<<endl;//DEBOGAGE
       }
      position=molpro_file.tellg();
      molpro_file.close();
       total+=sqrt(count[k]);
    //std::cout<<"probe2"<<std::endl;

       if(test2)
          break;
    }

    std::cout<<"AO S MATRIX READ IN MOLPRO OUTPUT WITHOUT ERROR"<<std::endl;
    position=0;
    
    MO_coeff_cation=new double[n_occ_tot**basis_size];
    for(int i=0;i!=n_occ_tot**basis_size;i++)
    {
       MO_coeff_neutral[i]=0;
       MO_coeff_cation[i]=0;
    }
    
   if(!search(&position, molpro_out_path,position, "NATURAL", 0, "ORBITALS"))
   {
      std::cout<<"LCAO COEFFICIENTS NOT FOUND IN MOLPRO OUTPUT"<<std::endl;
      return -1;
   }
       molpro_file.open(molpro_out_path.c_str());
       molpro_file.seekg(position);
           while(tmp_str != "Coefficients")
           {
              molpro_file>>tmp_str;
           }
           total=0;
           total2=0;
           for(int s=0;s!=n_sym;s++)
           {
            //cout<<"Coefficients of the neutral"<<endl;//DEBOGAGE
            for (int i=0; i!=2*sqrt(count[s]); i++)
            {
                molpro_file>>tmp_str;
    ;            //cout<<tmp_str<<endl;//DEBOGAGE
            }
            
            for (int j=total; j!=n_occ[s]+total; j++)
            {
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                for (int k=total2; k!=sqrt(count[s])+total2; k++)
                {
                    molpro_file>>floating;
                    MO_coeff_neutral[j**basis_size+k]=floating;
                    //cout<<MO_coeff_neutral[j**basis_size+k]<<"     ";//DEBOGAGE
                }//std::cout<<endl;//DEBOGAGE
            }//cout<<endl;
            total2+=sqrt(count[s]);
            total+=n_occ[s];
           }
           position=molpro_file.tellg();
           molpro_file.close();
    std::cout<<"MO COEFF NEUTRAL READ IN MOLPRO OUTPUT WITHOUT ERROR"<<std::endl;

   if(!search(&position, molpro_out_path,position, "NATURAL", 0, "ORBITALS"))
   {
      std::cout<<"LCAO COEFFICIENTS NOT FOUND IN MOLPRO OUTPUT"<<std::endl;
      return -1;
   }
       molpro_file.open(molpro_out_path.c_str());
       molpro_file.seekg(position);
           while(tmp_str != "Coefficients")
           {
              molpro_file>>tmp_str;
           }
           total=0;
           total2=0;
           for(int s=0;s!=n_sym;s++)
           {
            //cout<<"Coefficients of the cation"<<endl;//DEBOGAGE
            for (int i=0; i!=2*sqrt(count[s]); i++)
            {
                molpro_file>>tmp_str;
//                std::cout<<tmp_str<<std::endl;
            }
            for (int j=total; j!=n_occ[s]+total; j++)
            {
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                for (int k=total2; k!=sqrt(count[s])+total2; k++)
                {
                    molpro_file>>floating;
                    MO_coeff_cation[j**basis_size+k]=floating;
                }//cout<<endl;
            }
            total2+=sqrt(count[s]);
            total+=n_occ[s];
        }
    
        molpro_file.close();
    std::cout<<"MO COEFF CATION READ IN MOLPRO OUTPUT WITHOUT ERROR"<<std::endl;
   

   /*
        std::cout<<"#####################MO COeff NEUTRAL#####################"<<std::endl;
     for (int j=0; j!=n_occ_tot; j++)
     {
        std::cout<<"MO "<<j+1<<std::endl;
     for (int k=0; k!=*basis_size; k++)
     {
         cout<<MO_coeff_neutral[j**basis_size+k]<<"     ";//DEBOGAGE
     }std::cout<<endl;//DEBOGAGE
     }cout<<endl;
        std::cout<<"#########################################################"<<std::endl;
     //DEBOGAGE
    
        std::cout<<"#####################MO COeff cation#####################"<<std::endl;
     for (int j=0; j!=n_occ_tot; j++)
     {
        std::cout<<"MO "<<j+1<<std::endl;
     for (int k=0; k!=*basis_size; k++)
     {
         cout<<MO_coeff_cation[j**basis_size+k]<<"     ";//DEBOGAGE
     }std::cout<<endl;//DEBOGAGE
     }cout<<endl;
        std::cout<<"#########################################################"<<std::endl;
     //DEBOGAGE

        exit(EXIT_SUCCESS);
     */ 
    /*
    double answer(0);
    int mo_index(0);
    for(int i=0;i!=*basis_size;i++)
    {
       for(int j=0;j!=*basis_size;j++)
       {
          //answer+=MO_coeff_neutral[mo_index**basis_size+i]*MO_coeff_cation[mo_index**basis_size+j]*overlap[i**basis_size+j];
          answer+=MO_coeff_neutral[mo_index**basis_size+i]*MO_coeff_cation[mo_index**basis_size+j]*overlap[i**basis_size+j];
//          std::cout<<MO_coeff_neutral[mo_index**basis_size+i]<<"*"<<MO_coeff_cation[mo_index**basis_size+j]<<"*"<<overlap[i**basis_size+j]<<" = "<<MO_coeff_neutral[mo_index**basis_size+i]*MO_coeff_cation[mo_index**basis_size+j]*overlap[i**basis_size+j]<<" => "<<answer<<std::endl;
       }
    }
    std::cout<<" => "<<answer<<std::endl;
    exit(EXIT_SUCCESS);
    */
    temp=new double[*basis_size*n_occ_tot];
    temp2=new double[*basis_size*n_occ_tot];
    for(int i=0;i!=*basis_size*n_occ_tot;i++)
    {
       temp[i]=0;
       temp2[i]=0;
    }
    
    transpose(MO_coeff_neutral, temp2, n_occ_tot, *basis_size);
    

    matrix_product(temp, overlap, temp2, *basis_size, *basis_size, n_occ_tot);
    /*
    for (int i=0;i!=*basis_size; i++)
    {
        for (int j=0; j!=*n_occ; j++)
        {
            cout<<temp[i**n_occ+j]<<"    ";
        }cout<<endl;
    }cout<<endl;*///DEBOGAGE
    
    matrix_product(matrix2, MO_coeff_cation, temp, n_occ_tot, *basis_size, n_occ_tot);
    
    transpose(matrix2, matrix, n_occ_tot, n_occ_tot);
/*
    for(int i=0;i!=n_occ_tot;i++)
    {
       for(int j=0;j!=n_occ_tot;j++)
       {
          std::cout<<std::setw(15)<<std::setprecision(8)<<matrix[i*n_occ_tot+j];
       }std::cout<<std::endl<<std::endl;;
    }

    exit(EXIT_SUCCESS);
  */  
/*    
    double norm=0;
    for (int i=0; i!=*basis_size; i++)
    {
        for (int j=0; j!=*basis_size; j++)
        {
            norm+=MO_coeff_cation[i]*MO_coeff_cation[j]*overlap[i**basis_size+j];
        }
    }
    norm=sqrt(norm);
    cout<<"Norm of the MO is "<<norm<<endl<<endl;//DEBOGAGE
*/
    delete [] overlap;
    delete [] matrix2;
    delete [] MO_coeff_cation;
    delete [] temp;
    delete [] temp2;
    
    return 0;
}
int num_of_ci_reader(int *n_states_neut,int *n_states_cat,int *n_ci_neut,int *n_ci_cat,std::string file_address,int *n_occ,int n_sym=1)
{
    int position(0);
    std::string tmp_str("");
    int counter(0);
    double floating;
    int n_symocc(0);
    int test(0);
    bool test2(0);
    
    std::ifstream molpro_output;
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
   if(n_sym==1)
   {
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;
    
    while(tmp_str!="TOTAL")
    {
        molpro_output>>tmp_str;
        counter++;
    }
    
    *n_ci_neut=(counter-1)/(*n_states_neut+1);
    molpro_output.close();
    
    counter=0;
    
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;
    
    while(tmp_str!="TOTAL")
    {
        molpro_output>>tmp_str;
        counter++;
    }
    
    *n_ci_cat=(counter-1)/(*n_states_cat+1);
   }


   else
   {
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    n_symocc=0;
    for(int i=0;i!=n_sym;i++)
    {
       n_symocc+=bool(n_occ[i]);
    }
   //    std::cout<<"n_symocc= "<<n_symocc<<std::endl;

       int neut_states_sym=0;
       for(int i=0;i!=n_sym;i++)
       {
          neut_states_sym+=bool(n_states_neut[i]);
       }
    for(int i=0;i!=n_sym;i++)
    {
       test2=0;
       n_ci_neut[i]=0;
       if(n_states_neut[i]!=0)
       {
          counter=0;
          if(neut_states_sym!=1)
          {
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
          }
          molpro_output>>tmp_str;
    
          //std::cout<<"before neutral "<<tmp_str<<std::endl;//DEBOGAGE
          while(tmp_str!="TOTAL" && tmp_str!="CI" && tmp_str !="**********************************************************************************************************************************" )
          {
              molpro_output>>tmp_str;
            //  std::cout<<counter<<","<<tmp_str<<std::endl;
              counter++;
          }
          if(tmp_str == "TOTAL")
             test2=1;
          if(test2)
          {
             molpro_output>>tmp_str;
             for(int j=0;j!=n_states_neut[i];j++)
             {
                molpro_output>>tmp_str;
             }
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
          }
          else
             molpro_output>>tmp_str;
          //std::cout<<"after neutral "<<tmp_str<<std::endl;//DEBOGAGE
    
          n_ci_neut[i]=(counter-1)/(n_states_neut[i]+n_symocc);
       }
    }
    std::cout<<"PROBE 2 GETTING CI VECTOR SIZE FROM MOLPRO OUTPUT"<<std::endl;//DEBOGAGE
    
    position=molpro_output.tellg();
    molpro_output.close();
    
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    int cat_states_sym(0);
    for(int i=0;i!=n_sym;i++)
    {
       cat_states_sym+=bool(n_states_cat[i]);
    }
    for(int i=0;i!=n_sym;i++)
    {
       test2=0;
       n_ci_cat[i]=0;
       if(n_states_cat[i]!=0)
       {
          counter=0;

          if(cat_states_sym!=1)
          {
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
          }
          molpro_output>>tmp_str;
         
          //std::cout<<"PROBE LOOP CI VECTOR"<<std::endl;//DEBOGAGE

          //std::cout<<"before cation "<<tmp_str<<std::endl;//DEBOGAGE
          while(tmp_str!="TOTAL" && tmp_str!="CI" && tmp_str != "**********************************************************************************************************************************")
          {
            molpro_output>>tmp_str;
            counter++;
          }
          if(tmp_str == "TOTAL")
             test2=1;
          if(test2)
          {
              molpro_output>>tmp_str;
              //std::cout<<tmp_str<<std::endl;//DEBOGAGE
              for(int j=0;j!=n_states_cat[i];j++)
              {
                 molpro_output>>tmp_str;
                 //std::cout<<tmp_str<<std::endl;//DEBOGAGE
              }
              molpro_output>>tmp_str;
              // std::cout<<tmp_str<<std::endl;//DEBOGAGE
              molpro_output>>tmp_str;
              // std::cout<<tmp_str<<std::endl;//DEBOGAGE
              //std::cout<<"after cation "<<tmp_str<<std::endl;//DEBOGAGE
          } 
          else
             molpro_output>>tmp_str;
          //std::cout<<"counter is "<<counter<<" n_ci_cat is "<<(counter-1)/(n_states_cat[i]+n_symocc)<<std::endl;
          n_ci_cat[i]=(counter-1)/(n_states_cat[i]+n_symocc);

       }
    }
    }
    molpro_output.close();
    
    return 0;
}


int ci_vec_reader(int *n_states_neut_s,int *n_states_cat_s,int *n_occ,int *n_closed,int n_elec_neut,int *ci_size_neut,int *ci_size_cat,double **ci_vector_neut,double **ci_vector_cat,std::string file_address,int n_sym=1)
{
    int position(0);
    std::string tmp_str("");
    int elec_index(0);
    int state_index(0);
    double floating;
    int mo_index(0);
    int ci_index(0);
    int n_states_neut(0);
    int n_states_cat(0);
    int n_symocc(0);
    
    for(int i=0;i!=n_sym;i++)
    {
       n_states_neut+=n_states_neut_s[i];
       n_states_cat+=n_states_cat_s[i];
       std::cout<<n_states_cat_s[i]<<" cation states in sym "<<i<<std::endl;
    }
       std::cout<<n_states_cat<<" cation states "<<std::endl;

    std::ifstream molpro_output;

       std::cout<<"PROBE!!!!"<<std::endl;
    if(n_sym==1)
    {
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;

    for (int i=0; i!=*ci_size_neut; i++)
    {
        elec_index=0;
        molpro_output>>tmp_str;

        if(*n_closed!=0 && elec_index<2**n_closed)
        {
           for(int j=0;j!=*n_closed;j++)
           {
              ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
              ci_vector_neut[1][(n_elec_neut)*i+elec_index]=0;
              elec_index++;
              ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
              ci_vector_neut[1][(n_elec_neut)*i+elec_index]=1;
              elec_index++;

           }
        }
        for (int j=*n_closed; j!=*n_occ; j++)
        {
            if(tmp_str.at(j-*n_closed)=='0')
                continue;
            
            else if(tmp_str.at(j-*n_closed)=='2')
            {
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=0;
                elec_index++;
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=1;
            }
            
            else if(tmp_str.at(j-*n_closed)=='a')
            {
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=0;
            }
            
            else if(tmp_str.at(j-*n_closed)=='b')
            {
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=1;
            }
            
            elec_index++;
            
        }
        
        
        for (int j=0; j!=n_states_neut; j++)
        {
            molpro_output>>floating;
            ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+n_elec_neut+j]=floating;
            
        }
    }
    position=molpro_output.tellg();
    molpro_output.close();

    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;

    for (int i=0; i!=*ci_size_cat; i++)
    {
        elec_index=0;
        molpro_output>>tmp_str;

        if(*n_closed!=0 && elec_index<2**n_closed)
        {
           for(int j=0;j!=*n_closed;j++)
           {
              ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
              ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=0;
              elec_index++;
              ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
              ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=1;
              elec_index++;

           }
        }
        for (int j=*n_closed; j!=*n_occ; j++)
        {
            if(tmp_str.at(j-*n_closed)=='0')
                continue;
            
            else if(tmp_str.at(j-*n_closed)=='2')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=0;
                elec_index++;
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=1;
            }
            
            else if(tmp_str.at(j-*n_closed)=='a')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=0;
            }
            
            else if(tmp_str.at(j-*n_closed)=='b')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=1;
            }
            
            elec_index++;
            
        }
        
        
        for (int j=0; j!=n_states_cat; j++)
        {
            molpro_output>>floating;
            ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+n_elec_neut-1+j]=floating;
            
        }
    }
    }


    else
    {

       n_symocc=0;
       for(int i=0;i!=n_sym;i++)
       {
          n_symocc+=bool(n_occ[i]);
       }
       if(!search(&position, file_address,position, "CI", 0, "vector"))
       {
           std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
           return 1;
       }
        
       molpro_output.open(file_address.c_str());
    
       molpro_output.seekg(position);
    
       int neut_states_sym=0;
       for(int i=0;i!=n_sym;i++)
       {
          neut_states_sym+=bool(n_states_neut_s[i]);
       }
    
       for(int s=0;s!=n_sym;s++)
       {
         elec_index=0;
         if(n_states_neut_s[s]!=0)
         {
            if(neut_states_sym!=1)
            {
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
             molpro_output>>tmp_str;
            }
            molpro_output>>tmp_str;

            for (int i=0; i!=ci_size_neut[s]; i++)
            {
               mo_index=0;
               elec_index=0;
               for(int l=0;l!=n_symocc;l++)
               {
                  molpro_output>>tmp_str;
                  if(n_closed[l]!=0)
                  {
                     for(int j=0;j!=n_closed[l];j++)
                     {
                        ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                        ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=0;
                        elec_index++;
                        ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                        ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=1;
                        elec_index++;
                        mo_index++;
                     }
                 }
                 for(int j=0;j!=n_occ[l]-n_closed[l];j++)
                 {
                    if(tmp_str.at(j)=='0')
                    {
                       mo_index++;
                       continue;
                    }
            
                    else if(tmp_str.at(j)=='2')
                    {
                       ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                       ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=0;
                       elec_index++;
                       ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                       ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=1;
                       elec_index++;
                       mo_index++;
                       continue;
                    }
            
                    else if(tmp_str.at(j)=='a' || tmp_str.at(j)=='/')
                    {
                       ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                       ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=0;
                       elec_index++;
                       mo_index++;
                       continue;
                    }
             
                    else if( tmp_str.at(j) == 'b' || tmp_str.at(j) == '\\')
                    {
                       ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                       ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=1;
                       elec_index++;
                       mo_index++;
                       continue;
                    }
                  else 
                  {
                     std::cout<<"ERROR unrecognized symbol "<<tmp_str.at(j)<<std::endl;
                     std::cout<<"EXIT"<<std::endl;
                     exit(EXIT_SUCCESS);
                  }
                }
            }
            elec_index=0;
            state_index=0;
            for(int w=0;w!=s;w++)
            {
               state_index+=n_states_neut_s[w];
            }
        
            for (int j=state_index; j!=n_states_neut_s[s]+state_index; j++)
            {
                molpro_output>>floating;
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+n_elec_neut+j]=floating;
            
            }
            ci_index++;
        }
      
        molpro_output>>tmp_str;
        if(tmp_str == "TOTAL")
        {
           molpro_output>>tmp_str;
           for(int w=0;w!=n_states_neut_s[s];w++)
           {
              molpro_output>>tmp_str;
           }
           molpro_output>>tmp_str;
           molpro_output>>tmp_str;
        }
        else
           molpro_output>>tmp_str;
    }
    }
    position=molpro_output.tellg();
    molpro_output.close();
    
   ci_index=0;
   elec_index=0;
   state_index=0;
    
   std::cout<<"Neutral part of CI vector read. Now reading cation CI vector"<<std::endl;

    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);

    int cat_states_sym(0);
    for(int i=0;i!=n_sym;i++)
    {
       cat_states_sym+=bool(n_states_cat_s[i]);
    }
    std::cout<<cat_states_sym<<" symmetry in cationic states"<<std::endl;//DEBOGAGE
    
    for(int s=0;s!=n_sym;s++)
    {
       //std::cout<<"probe sym "<<s<<std::endl;//DEBOGAGE
      if(n_states_cat_s[s]!=0)
      {
        //std::cout<<"probe loop sym"<<std::endl;//DEBOGAGE
       if(cat_states_sym!=1) //!!!!!!! THIS SHOULD BE ONLY COMMENTED IF WE WANT A SINGLE STATE OF THE CAITON BUT THE MOLPRO COMPUTATION WAS DONE ON SEVERAL STATES
       {
          molpro_output>>tmp_str;
          molpro_output>>tmp_str;
          molpro_output>>tmp_str;
          molpro_output>>tmp_str;
       }
        molpro_output>>tmp_str;
    
        for (int i=0; i!=ci_size_cat[s]; i++)
        {
           elec_index=0;
           mo_index=0;
           for(int l=0;l!=n_symocc;l++)
           {
              molpro_output>>tmp_str;
               if(n_closed[l]!=0)
               {
                   for(int j=0;j!=n_closed[l];j++)
                   {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=0;
                      elec_index++;
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=1;
                      elec_index++;
                      mo_index++;
                   }
               }
               for(int j=0;j!=n_occ[l]-n_closed[l];j++)
               {
                  if(tmp_str.at(j)=='0')
                  {
                      mo_index++;
                  //    continue;
                  }
            
                  else if(tmp_str.at(j)=='2')
                  {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=0;
                      elec_index++;
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=1;
                      elec_index++;
                      mo_index++;
                  //    continue;
                  }
            
                  else if(tmp_str.at(j)=='a' || tmp_str.at(j)=='/')
                  {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=0;
                      elec_index++;
                      mo_index++;
                  //    continue;
                  }
            
                  else if( tmp_str.at(j) == 'b' || tmp_str.at(j) == '\\')
                  {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=1;
                      elec_index++;
                      mo_index++;
                  //    continue;
                  }
                  else 
                  {
                     std::cout<<"ERROR unrecognized symbol "<<tmp_str.at(j)<<std::endl;
                     std::cout<<"EXIT"<<std::endl;
                     exit(EXIT_SUCCESS);
                  }

               }
         }
         /*std::cout<<"Configuration "<<ci_index<<" state sym "<<s<<std::endl;
         for(int j=0;j!=n_elec_neut-1;j++)
         {
            std::cout<<int(ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+j])<<" ";
         }std::cout<<std::endl;*/
           //std::cout<<"probe 2"<<std::endl;//DEBOGAGE

           elec_index=0;
           state_index=0;
         for(int w=0;w!=s;w++)
         {
              state_index+=n_states_cat_s[w];
         }
         for (int j=state_index; j!=n_states_cat_s[s]+state_index; j++)
         {
                molpro_output>>floating;
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+n_elec_neut-1+j]=floating;
         }
         ci_index++;
        }
        molpro_output>>tmp_str;
        if(tmp_str == "TOTAL")
        {
           molpro_output>>tmp_str;
           for(int w=0;w!=n_states_cat_s[s];w++)
           {
              molpro_output>>tmp_str;
           }
           molpro_output>>tmp_str;
           molpro_output>>tmp_str;
        }
        else
           molpro_output>>tmp_str;
        //std::cout<<"probe loop sym"<<std::endl;//DEBOGAGE
   }
    }
}
          /// std::cout<<"probe end"<<std::endl;//DEBOGAGE
    
    molpro_output.close();
    return 0;
}

bool basis_size_data_reader(int n_sym, int* basis_size_sym,int** contraction_number,std::string file_address)
{
   using namespace std;
   ifstream input;
   int count(0);
   stringstream sstream;
   string teststring;
   int basis_size(0);
   for(int i=0;i!=n_sym;i++)
   {
      basis_size+=basis_size_sym[i];
   }
   string temps;
   int pos(0);
   if(!search(&pos,file_address.c_str(),0,"BASIS",0,"DATA"))
   {
      std::cout<<"BASIS DATA NOT FOUND IN MOLPRO OUTPUT FILE"<<std::endl;
      exit(EXIT_FAILURE);
   }

   input.open(file_address.c_str());

   if(!input.is_open())
   {
      std::cout<<"ERROR WHILE OPENING INPUT FILE "<<std::endl<<file_address.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      input.seekg(pos);
      for(int i=0;i!=7;i++)
      {
         input>>temps;
      }
      for(int s=0;s!=n_sym;s++)
      {
         //std::cout<<basis_size_sym[s]<<std::endl;
         for(int i=0;i!=basis_size_sym[s];i++)
         {
            for(int j=0;j!=4;j++)
            {
               input>>temps; 
//               std::cout<<temps<<std::endl;
            }
            temps="";
            sstream.str("");
            if(i<basis_size_sym[s]-1)
            {
//               std::cout<<"!"<<i<<std::endl;
               sstream<<i+2<<"."<<s+1;
            }
            else if(s<n_sym-1)
            {
               sstream<<"1."<<s+2;
            }
            else if (i==basis_size_sym[s]-1 && s == n_sym-1)
            {
               sstream<<i+1<<"."<<n_sym;
            }
            teststring=sstream.str();

//            std::cout<<teststring<<"±±±±±±±±±±"<<std::endl;
            count=0;
            while(temps!=teststring && temps != "NUCLEAR")
            {
               input>>temps;
               count++;
//               std::cout<<temps.c_str()<<"!="<<teststring.c_str()<<"!!!"<<bool(temps!=teststring)<<std::endl;
               if(input.eof())
               {
                  break;
               }
            }
               if(input.eof())
               {
                  std::cout<<"end of loop"<<std::endl;
                  break;
               }
//            std::cout<<temps<<std::endl;
            contraction_number[s][i]=count/2;
//            std::cout<<count<<"=>"<<contraction_number[s][i]<<std::endl;

         }
      }
      input.close();
}
return 0;
}
bool basis_data_reader(int n_sym, int* basis_size_sym,int** contraction_number,double*** contraction_coeff,double*** contraction_zeta,int** nucl_basis_func,std::string **basis_func_type,std::string file_address,int ***angular_mom_numbers)
{
   using namespace std;
   ifstream input;
   int count(0);
   stringstream sstream;
   string teststring;
   int basis_size(0);
   double temp_norm(0);
   bool test(0);
   int l(0);
   int l_val(0);
   for(int i=0;i!=n_sym;i++)
   {
      basis_size+=basis_size_sym[i];
   }
   string temps;
   int pos(0);
   search(&pos,file_address.c_str(),0,"BASIS",0,"DATA");

   input.open(file_address.c_str());

   if(!input.is_open())
   {
      std::cout<<"ERROR WHILE OPENING INPUT FILE "<<std::endl<<file_address.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      input.seekg(pos);
      for(int i=0;i!=7;i++)
      {
         input>>temps;
      }
      for(int s=0;s!=n_sym;s++)
      {
         for(int i=0;i!=basis_size_sym[s];i++)
         {
            input>>temps;
//            std::cout<<"1  "<<temps<<std::endl;
            input>>temps;
            input>>temps;
            nucl_basis_func[s][i]=atoi(temps.c_str());
            std::cout<<i+1<<"."<<s+1<<"- nucleus "<<nucl_basis_func[s][i]<<std::endl;
            input>>temps;
            basis_func_type[s][i]=temps.c_str();
            angular_mom_numbers[s][i][0]=l_number(basis_func_type[s][i]);
            angular_mom_numbers[s][i][1]=ml_number(basis_func_type[s][i],angular_mom_numbers[s][i][0]);
//            std::cout<< basis_func_type[s][i].c_str()<<" -> l = "<<angular_mom_numbers[s][i][0]<<" and ml = "<<angular_mom_numbers[s][i][1]<<std::endl;
            l_val=angular_mom_numbers[s][i][0];
            test=0;
            l=0;
            while(test==0)
            {
//               std::cout<<l<<"/"<<contraction_number[s][i]<<std::endl;
               input>>temps;
               contraction_zeta[s][i][l]=atof(temps.c_str());
               input>>temps;
               contraction_coeff[s][i][l]=atof(temps.c_str());
               contraction_coeff[s][i][l]/=(sqrt(0.5*intplushalf_gamma(l_val+1)/pow(2*contraction_zeta[s][i][l],l_val+1.5)));
               l++;
               if(l>=contraction_number[s][i])
               {
                  test=1;
               }

            }
//            std::cout<<"2  "<<temps<<std::endl;
         }
      }
      input.close();
   }
 /*  for(int s=0;s!=n_sym;s++)
   {
      for(int i=0;i!=basis_size_sym[s];i++)
      {
         l_val=l_number(basis_func_type[s][i].c_str());
         temp_norm=0;
         for(l=0;l!=contraction_number[s][i];l++)
         {
         //   std::cout<<contraction_coeff[s][i][l]<<"!!!"<<l<<"$$$"<<i<<"###"<<s<<std::endl;
            for(int lp=0;lp!=contraction_number[s][i];lp++)
            {
               temp_norm+=0.5*contraction_coeff[s][i][l]*contraction_coeff[s][i][lp]*intplushalf_gamma(l_val+1)/pow(contraction_zeta[s][i][l]+contraction_zeta[s][i][lp],l_val+1.5);
//               std::cout<<lp<<"  "<<l_val<<"   "<<intplushalf_gamma(l_val)<<"   "<<std::endl;
            }
         }
            std::cout<<"norm of basis function "<<i<<" symmetry "<<s<<":"<<temp_norm<<std::endl;

         if(temp_norm==0)
         {
            std::cout<<"ERROR ! Norm of basis function "<<i<<" = 0 "<<std::endl;
            exit(EXIT_FAILURE);
         }
         for(l=0;l!=contraction_number[s][i];l++)
         {
            contraction_coeff[s][i][l]/=sqrt(temp_norm);
            //std::cout<<"==>"<<contraction_coeff[s][i][l]<<"!!!"<<l<<"$$$"<<i<<"###"<<s<<std::endl;
         }
         
         temp_norm=0;
         for(l=0;l!=contraction_number[s][i];l++)
         {
            for(int lp=0;lp!=contraction_number[s][i];lp++)
            {
               temp_norm+=contraction_coeff[s][i][l]*contraction_coeff[s][i][lp]*pow(acos(-1)/(contraction_zeta[s][i][l]+contraction_zeta[s][i][lp]),1.5);
            }
         }
            std::cout<<"New norm of basis function "<<i<<" symmetry "<<s<<":"<<temp_norm<<std::endl;
            
      }
   }*/
      return 0;

}

bool search(int *match_loc,std::string file_address,int research_from,std::string pattern1, int num_of_entry_between_patterns12,std::string pattern2,int num_of_entry_between_patterns23,std::string pattern3)
{
    using namespace std;
    string tmp_str;
    ifstream molpro_output;
    int count(0);
    bool test(0);
    
    molpro_output.open(file_address.c_str());
    if (!molpro_output.is_open())
    {
        std::cout<<"Issue while opening molpro output file"<<endl;
        return 0;
    }
    else
    {
        molpro_output.seekg(research_from);
        
        molpro_output>>tmp_str;
        do
        {
            while (tmp_str!=pattern1)
            {
                
                molpro_output>>tmp_str;
                
                if(molpro_output.eof())
                {
                    molpro_output.close();
                    return 0;
                }
               // std::cout<<tmp_str<<" (1) = "<<pattern1<<std::endl;//DEBOGAGE
            }
            
            if(pattern2!="")
            {
                
                
                do
                {
                    molpro_output>>tmp_str;
                    count++;
                }while(count<=num_of_entry_between_patterns12);
                count=0;
                
                //std::cout<<tmp_str<<" (2) = "<<pattern2<<std::endl;//DEBOGAGE
                
                if (tmp_str==pattern2)
                {
                    if(pattern3=="")
                    {       
                //       std::cout<<"!probe!!"<<std::endl;
                        *match_loc=int(molpro_output.tellg());
                        molpro_output.close();
                        test=1;
                        return 1;
                    }
                    else
                    {
                       do
                       {
                          molpro_output>>tmp_str;
                          count++;
                       }while(count<=num_of_entry_between_patterns23);
                       count=0;

                       if (tmp_str==pattern3)
                       {
                          //std::cout<<"position = "<<filestream.pos()<<std::endl;//DEBOGAGE
                          *match_loc=int(molpro_output.tellg());
                          molpro_output.close();
                          test=1;
                          return 1;
                       }
                    }
                }

            }
            else
            {
               *match_loc=int(molpro_output.tellg());
               molpro_output.close();
               test=1;
               return 1;
            }
        }while(test!=1);//Search for the entry in the file

    }
    molpro_output.close();
    return 1;
}


charTypeT charType(char c)
{
    if(isdigit(c))return digit;
    if(isalpha(c))return alpha;
    return other;
}

std::string separateThem(std::string inString)
{
   using namespace std;
  string oString = "";charTypeT st=other;
    for(int i=0;i!=inString.size();i++)
    {
       char c = inString[i];
        if( (st==alpha && charType(c)==digit) || (st==digit && charType(c)==alpha) )
          return oString;
         // oString.push_back(' ');
        oString.push_back(c);st=charType(c);
    }
}
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
