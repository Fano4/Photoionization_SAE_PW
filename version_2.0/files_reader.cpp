
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
            return -1;
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

int overlap_MO(double matrix[],int* n_occ,int* basis_size,std::string molpro_out_path,int n_sym=1)
{
    bool test(0);
    bool test2(0);
    int *count=new int [n_sym];
    int total(0);
    int total2(0);
    int position(0);
    double *overlap;
    double *matrix2;
    double *MO_coeff_neutral;
    double *MO_coeff_cation;
    double floating;
    double *temp;
    double *temp2;
    int n_occ_tot(0);

    for(int k=0;k!=n_sym;k++)
    {
       n_occ_tot+=n_occ[k];
    } 
    std::cout<<n_occ_tot<<std::endl;
    using namespace std;
    
    string tmp_str;
    
    
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

       for (int i=total; i<sqrt(count[k])+total; i++)
       {
          for(int j=total;j<sqrt(count[k])+total;j++)
          {
             molpro_file>>floating;
             overlap[i**basis_size+j]=floating;
          }
                     //cout<<overlap[i**basis_size+total]<<endl;//DEBOGAGE
       }
       position=molpro_file.tellg();
       molpro_file.close();
       total+=sqrt(count[k]);

       if(test2)
          break;
    }

    position=0;
    
    MO_coeff_neutral=new double[n_occ_tot**basis_size];
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
           for(int i=0;i!=7;i++)
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
                //cout<<tmp_str<<endl;//DEBOGAGE
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

   if(!search(&position, molpro_out_path,position, "NATURAL", 0, "ORBITALS"))
   {
      std::cout<<"LCAO COEFFICIENTS NOT FOUND IN MOLPRO OUTPUT"<<std::endl;
      return -1;
   }
       molpro_file.open(molpro_out_path.c_str());
       molpro_file.seekg(position);
           for(int i=0;i!=7;i++)
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


    
    
   /* double norm=0;
    for (int i=0; i!=*basis_size; i++)
    {
        for (int j=0; j!=*basis_size; j++)
        {
            norm+=MO_coeff_cation[i]*MO_coeff_cation[j]*overlap[i**basis_size+j];
        }
    }
    norm=sqrt(norm);
    cout<<"Norm of the MO is "<<norm<<endl<<endl;*///DEBOGAGE

    delete [] overlap;
    delete [] matrix2;
    delete [] MO_coeff_neutral;
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
       std::cout<<"n_symocc= "<<n_symocc<<std::endl;

       int neut_states_sym=0;
       for(int i=0;i!=n_sym;i++)
       {
          neut_states_sym+=bool(n_states_neut[i]);
       }
    for(int i=0;i!=n_sym;i++)
    {
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
       while(tmp_str!="TOTAL")
       {
           molpro_output>>tmp_str;
           counter++;
       }
        molpro_output>>tmp_str;
        for(int j=0;j!=n_states_neut[i];j++)
        {
           molpro_output>>tmp_str;
        }
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
        //std::cout<<"after neutral "<<tmp_str<<std::endl;//DEBOGAGE
    
       n_ci_neut[i]=(counter-1)/(n_states_neut[i]+n_symocc);
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
    
    int cat_states_sym(0);
    for(int i=0;i!=n_sym;i++)
    {
       cat_states_sym+=bool(n_states_cat[i]);
    }
    for(int i=0;i!=n_sym;i++)
    {
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

       //std::cout<<"before cation "<<tmp_str<<std::endl;//DEBOGAGE
        while(tmp_str!="TOTAL")
        {
            molpro_output>>tmp_str;
            counter++;
        }
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
    }
    std::ifstream molpro_output;

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
        
        elec_index=0;
        
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
        
        elec_index=0;
        
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
            
                  else if(tmp_str.at(j)=='a')
                  {
                      ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                      ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=0;
                      elec_index++;
                      mo_index++;
                      continue;
                   }
            
                  else if(tmp_str.at(j)=='b')
                  {
                      ci_vector_neut[0][(n_elec_neut+n_states_neut)*ci_index+elec_index]=mo_index;
                      ci_vector_neut[1][(n_elec_neut)*ci_index+elec_index]=1;
                      elec_index++;
                      mo_index++;
                      continue;
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
        molpro_output>>tmp_str;
        for(int w=0;w!=n_states_neut_s[s];w++)
        {
           molpro_output>>tmp_str;
        }
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
    }
    }
    position=molpro_output.tellg();
    molpro_output.close();
    
   ci_index=0;
    
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
    
    for(int s=0;s!=n_sym;s++)
    {
      if(n_states_cat_s[s]!=0)
      {
       if(cat_states_sym!=1)
       {
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
       }
        molpro_output>>tmp_str;
    
        for (int i=0; i!=ci_size_cat[s]; i++)
        {
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
                      continue;
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
                      continue;
                  }
            
                  else if(tmp_str.at(j)=='a')
                  {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=0;
                      elec_index++;
                      mo_index++;
                      continue;
                  }
            
                  else if(tmp_str.at(j)=='b')
                  {
                      ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*ci_index+elec_index]=mo_index;
                      ci_vector_cat[1][(n_elec_neut-1)*ci_index+elec_index]=1;
                      elec_index++;
                      mo_index++;
                      continue;
                  }

               }
         }
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
        molpro_output>>tmp_str;
        for(int w=0;w!=n_states_cat_s[s];w++)
        {
           molpro_output>>tmp_str;
        }
        molpro_output>>tmp_str;
        molpro_output>>tmp_str;
   }
    }
}
    
    molpro_output.close();
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
                //std::cout<<tmp_str<<" (1) = "<<pattern1<<std::endl;//DEBOGAGE
            }
            
            if(pattern2!="")
            {
                
                
                do
                {
                    molpro_output>>tmp_str;
                    count++;
                }while(count<=num_of_entry_between_patterns12);
                count=0;
                
                //std::cout<<tmp_str.toStdString()<<" (2) = "<<pattern2<<std::endl;//DEBOGAGE
                
                if (tmp_str==pattern2)
                {
                    if(pattern3=="")
                    {
                        
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


int potential_and_dipole_reader(std::string file_address,int ini_pos,int n_states, double *pot, double *dipole_x, double *dipole_y, double *dipole_z)
{
   using namespace std;
   int found_pos(ini_pos);

   std::stringstream file_indenter;
   std::string variable_pattern;

   for(int i=0;i!=n_states;i++)
   {
      file_indenter.str("");
      file_indenter<<i+1<<"."<<1;
      variable_pattern=file_indenter.str();

      if(!search(&found_pos,file_address,ini_pos,"!MCSCF",1,variable_pattern,0,"Energy"))
      {
         std::cout<<" POTENTIAL NOT FOUND IN MOLPRO OUTPUT"<<std::endl<<"STATE "<<variable_pattern<<std::endl;
         return 0;
      }
      else
      {
         ifstream  molpro_output;
         molpro_output.open(file_address.c_str());
         molpro_output.seekg(found_pos);
         molpro_output>>pot[i];
         molpro_output.close();
      }


      for(int j=i;j!=n_states;j++)
      {

         file_indenter.str("");
         file_indenter<<"<"<<j+1<<"."<<1<<"|DMX|"<<i+1<<".1>";
         variable_pattern=file_indenter.str();

         if(!search(&found_pos,file_address,ini_pos,"!MCSCF",1,variable_pattern))
         {
            dipole_x[(n_states * i) + j - ((i * (i+1)) / 2)]=0;
         }
         else
         {
            ifstream  molpro_output;
            molpro_output.open(file_address.c_str());
            molpro_output.seekg(found_pos);
            molpro_output>>dipole_x[(n_states * i) + j - ((i * (i+1)) / 2)];
            molpro_output.close();
         }

         file_indenter.str("");
         file_indenter<<"<"<<j+1<<"."<<1<<"|DMY|"<<i+1<<".1>";
         variable_pattern=file_indenter.str();

         if(!search(&found_pos,file_address,ini_pos,"!MCSCF",1,variable_pattern))
         {
            dipole_y[(n_states * i) + j - ((i * (i+1)) / 2)]=0;
         }
         else
         {
            ifstream  molpro_output;
            molpro_output.open(file_address.c_str());
            molpro_output.seekg(found_pos); 
            molpro_output>>dipole_y[(n_states * i) + j - ((i * (i+1)) / 2)];
            molpro_output.close();
         }

         file_indenter.str("");
         file_indenter<<"<"<<j+1<<"."<<1<<"|DMZ|"<<i+1<<".1>";
         variable_pattern=file_indenter.str();

         if(!search(&found_pos,file_address,ini_pos,"!MCSCF",1,variable_pattern))
         {
            dipole_z[(n_states * i) + j - ((i * (i+1)) / 2)]=0;
         }
         else
         {
            ifstream  molpro_output;
            molpro_output.open(file_address.c_str());
            molpro_output.seekg(found_pos); 
            molpro_output>>dipole_z[(n_states * i) + j - ((i * (i+1)) / 2)];
            molpro_output.close();
         }

      }
   }
   double pot_gs=pot[0];

   for(int i=0;i!=n_states;i++)
   {
      pot[i]-=pot_gs;
      std::cout<<pot[i]<<std::endl;
   }


   return found_pos;
}


bool efield_param_reader(std::string file_address,int *n_pulses,double **energy,double **intensity,double **origin, double **sigma,double **CEP,int **ppol)
{
   using namespace std;
   int match_loc(0);
   string pattern;

   pattern="$$n_pulses";   

   if(!search(&match_loc,file_address,match_loc,pattern))
   {
      cout<<"NUMBER OF PULSES NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl;
      return 0;
   }
   else
   {
      ifstream efield_input;
      efield_input.open(file_address.c_str());
      efield_input>>*n_pulses;
      efield_input.close();
   }

   *energy=new double[*n_pulses];
   *intensity=new double[*n_pulses]; 
   *origin=new double[*n_pulses]; 
   *sigma=new double[*n_pulses]; 
   *CEP=new double[*n_pulses]; 

   pattern="$$energy";

   match_loc=0;
   for(int i=0;i!=*n_pulses;i++)
   {
     if(!search(&match_loc,file_address,match_loc,pattern))
     {   
        cout<<"PULSE ENERGY PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl;
        return 0;
     }    
     else
     {   
        ifstream efield_input;
        efield_input.open(file_address.c_str());
        efield_input>>*energy[i];
        efield_input.close();
     }   

    }

   pattern="$$intensity";                                                         

   match_loc=0;                                                                
   for(int i=0;i!=*n_pulses;i++)                                               
   {
     if(!search(&match_loc,file_address,match_loc,pattern))          
     {                                                                         
        cout<<"PULSE INTENSITY PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl;                                               
        return 0;                                                              
     }                                                                         
     else                                                                      
     {                                                                         
        ifstream efield_input;                                                
        efield_input.open(file_address.c_str());           
        efield_input>>*intensity[i];                                             
        efield_input.close();                                                 
     }   

    }

    pattern="$$origin";
   
   match_loc=0;
   for(int i=0;i!=*n_pulses;i++)                                               
   { 
     if(!search(&match_loc,file_address,match_loc,pattern))          
     {
        cout<<"PULSE ORIGIN PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl<<"PULSE "<<i<<std::endl;
        return 0;
     }
     else
     {
        ifstream efield_input;
        efield_input.open(file_address.c_str());
        efield_input>>*origin[i];
        efield_input.close();
     }

    }

    pattern="$$sigma";
   
   match_loc=0;
   for(int i=0;i!=*n_pulses;i++)
   {
     if(!search(&match_loc,file_address,match_loc,pattern))
     {
        cout<<"PULSE SIGMA PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl<<"PULSE "<<i<<std::endl;
        return 0;
     } 
     else
     {
        ifstream efield_input;
        efield_input.open(file_address.c_str());
        efield_input>>*sigma[i];
        efield_input.close();
     }  
     
    }
   pattern="$$CEP";
   
   match_loc=0;
   for(int i=0;i!=*n_pulses;i++)
   {   
     if(!search(&match_loc,file_address,match_loc,pattern))
     {   
        cout<<"PULSE CEP PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl<<"PULSE "<<i<<std::endl;
        return 0;
     }   
     else
     {   
        ifstream efield_input;
        efield_input.open(file_address.c_str());
        efield_input>>*CEP[i];
        efield_input.close();
     }   
    
    } 
    pattern="$$polarisation";

   match_loc=0;
   for(int i=0;i!=*n_pulses;i++)
   {
     if(!search(&match_loc,file_address,match_loc,pattern))
     {
        cout<<"PULSE POLARISATION PARAMETER NOT FOUND IN ELECTRIC FIELD INPUT :"<<std::endl<<file_address<<std::endl<<"PULSE "<<i<<std::endl;
        return 0;
     }
     else
     {
        ifstream efield_input;
        efield_input.open(file_address.c_str());
        efield_input>>*ppol[i];
        efield_input.close();
     }

    }
   return 1;
}
int photoion_input_reader(std::string input_loc, elec_struct *elec_target,cube_struct *cube_target, continuum_struct *continuum_target)
{
   using namespace std;
   ifstream input;
   input.open(input_loc.c_str());

   //electronic structure variables 
   bool molpro_input_loc_checked(0);
   bool mo_cube_loc_set(0);
   string mo_cube_loc;

   //cube_struct_variables
   bool cube_dim_set(0);
   bool nx_set(0);
   bool ny_set(0);
   bool nz_set(0);
   bool xmin_set(0);
   bool xmax_set(0);
   bool ymin_set(0);
   bool ymax_set(0);
   bool zmin_set(0);
   bool zmax_set(0);
   bool cube_titles_set(0);
   bool cube_geom_set(0);
   bool num_nucl_set(0);
   int nx,ny,nz(0);
   double xmin,ymin,zmin(0);
   double xmax,ymax,zmax(0);
   int num_nucl(0);
   int *nucl_Z(NULL);
   double **nucl_pos(NULL);

   //continuum struct variables
   bool k_range_set(0);
   bool grid_set(0);
   bool kmin_set(0);
   bool kmax_set(0);
   bool nk_set(0);
   bool ntheta_set(0);
   bool nphi_set(0);
   bool random_grid(0);
   bool random_grid_set(0);

   double kmin,kmax(0);
   int nk,ntheta,nphi(0);
   int n_points_sphere(0);

   //utility variables
   char c_comment;
   string l_comment;
   int position(0);
   string parameter_label;
   int parameter_int;
   double parameter_double;
   string parameter_string;
   int num_of_atom(0);

   if(!input.is_open())
   {
      cout<<"INPUT FILE DOES NOT EXIST"<<endl<<input_loc.c_str()<<endl;
      return 0;
   }
   else
   { 
      do
      {
         c_comment=input.get();
         if(c_comment=='!')
         {
            //ignore comment lines
            getline(input,l_comment);
            continue;
         }
         else if(c_comment==' ')
         {
            //ignore void characters at the beginning of a line
            continue;
         }
         else
         {
            //take one step back
            position=input.tellg();
            input.seekg(position-1);
            //read the parameter label. The label is identified using a switch condition
            input>>parameter_label; 
                  //electronic structure input section
                  if(parameter_label == "molpro_file")
                  {
                     input>>input_loc;
                  }
                  else if(parameter_label == "MO_cubes_loc")
                  {
                     input>>parameter_string;
                     elec_target->set_mo_cube(parameter_string); 
                  }
                  //cube structure section
                  else if(parameter_label == "grid_size_x")
                  {
                     input>>nx; 
                     if(nx<=0)
                     {
                        cout<<"ERROR: NEGATIVE OR ZERO GRID SIZE MAKES NO SENSE. PLEASE CHECK INPUT"<<endl;
                        return 8;
                     }
                     nx_set=1;
                  }
                  else if(parameter_label == "grid_size_y")
                  {
                     input>>ny;
                     if(ny<=0)
                     {
                        cout<<"ERROR: NEGATIVE OR ZERO GRID SIZE MAKES NO SENSE. PLEASE CHECK INPUT"<<endl;
                        return 8;
                     }
                     ny_set=1;
                  }
                  else if(parameter_label == "grid_size_z")
                  {
                     input>>nz;
                     nz_set=1;
                     if(nz<=0)
                     {
                        cout<<"ERROR: NEGATIVE OR ZERO GRID SIZE MAKES NO SENSE. PLEASE CHECK INPUT"<<endl;
                        return 8;
                     }
                  }
                  else if(parameter_label == "min_x")
                  {
                     input>>xmin;
                     xmin_set=1;
                  }
                  else if(parameter_label == "min_y")
                  {
                     input>>ymin;
                     ymin_set=1;
                  }
                  else if(parameter_label == "min_z")
                  {
                     input>>zmin;
                     zmin_set=1;
                  }
                  else if(parameter_label == "max_x")
                  {
                     input>>xmax;
                     xmax_set=1;
                  }
                  else if(parameter_label == "max_y")
                  {
                     input>>ymax;
                     ymax_set=1;
                  }
                  else if(parameter_label == "max_z")
                  {
                     input>>zmax;
                     zmax_set=1;
                  }
                  //molecular Geometry section
                  else if(parameter_label == "num_of_nucl")
                  {
                     input>>num_nucl;
                     nucl_Z=new int[num_nucl];
                     nucl_pos=new double*[num_nucl];
                     for(int i=0;i!=num_nucl;i++)
                     {
                        nucl_pos[i]=new double[3];
                     }
                     num_nucl_set=1;
                  }
                  else if(parameter_label == "atom")
                  {
                     if(!num_nucl_set)
                     {
                        cout<<"ERROR: NUMBER OF NUCLEI NOT SET BEFORE GEOMETRY SPECIFICATION"<<endl;
                        return 8;
                     }
                     if(num_of_atom+1>num_nucl)
                     {
                        cout<<"ERROR: TOO MANY ATOMS IN GEOMETRY INPUT. EXPECTED "<<num_nucl<<", DETECTED "<<num_of_atom+1<<". PLEASE CHECK INPUT"<<endl;
                        return 8;
                     }
                     input>>nucl_Z[num_of_atom];
                     input>>nucl_pos[num_of_atom][0];
                     input>>nucl_pos[num_of_atom][1];
                     input>>nucl_pos[num_of_atom][2];
                     num_of_atom++;
                  }
                  //Continuum subspace section
                  else if(parameter_label == "kmin")
                  {
                     input>>kmin;
                     kmin_set=1;
                  }
                  else if(parameter_label == "kmax")
                  {
                     input>>kmax;
                     kmax_set=1;
                  }
                  else if(parameter_label == "grid_size_k")
                  {
                     input>>nk;
                     nk_set=1;
                  }
                  else if(parameter_label == "random_grid")
                  {
                     input>>random_grid;
                  }
                  else if(parameter_label == "n_points_sphere")
                  {
                     input>>n_points_sphere;
                     random_grid_set=1;
                  }
                  else if(parameter_label == "grid_size_polar")
                  {
                     input>>ntheta;
                     ntheta_set=1;
                  }
                  else if(parameter_label == "grid_size_azim")
                  {
                     input>>nphi;
                     nphi_set=1;
                  }
                  else
                  {
                     std::cout<<"PARAMETER NOT RECOGNIZED IN INPUT FILE :"<<parameter_label.c_str()<<std::endl;
                     return 9;
                  }
            //Check parameters setting
            if(nx_set && ny_set && nz_set && xmin_set && xmax_set && ymin_set && ymax_set && zmin_set && zmax_set)
               cube_dim_set=1;
            if(num_of_atom==num_nucl)
               cube_geom_set=1;
            if(kmin_set && kmax_set)
               k_range_set=1;
            if(nk_set && ((random_grid && random_grid_set) || (ntheta_set && nphi_set)))
               grid_set=1;

         }
      }while(!input.eof());
      input.close();
//Set up the target objects
      if(cube_dim_set)
         cube_target->set_cube_dim(nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax);
   }
   return 1;
}
