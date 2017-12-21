//
//  molpro_out_reader.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 19/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#include "molpro_out_reader.hpp"

// This block of functions reads molpro output and needs the following gprints: basis,overlap S,CI-vectors,orbitals (occ)

//It begins by computing the determinants of the overlap between the cationic and the neutral molecular orbitals for every possible ionization from neutral.

int size_query(int* n_occ,int* basis_size,std::string molpro_out_path)
{
    bool test(0);
    
    using namespace std;
    
    string tmp_str;
    int integer;
    double floating;
    
    
    
        ifstream molpro_file;
        molpro_file.open(molpro_out_path.c_str());
        
        if (!molpro_file.is_open())
        {
            cout<<"ERROR DURING MOLPRO OUTPUT FILE ACQUISITION"<<endl<<molpro_out_path.c_str()<<endl;
            return -1;
        }
        
        
        
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
                
                //cout<<"Basis set size found!! Size= "<<basis_size<<endl;//DEBOGAGE
                test=1;
            }
        }while(test!=1);//Search for the basis set size in molpro output
    
    do
    {
        while (tmp_str!="Number")
        {
            molpro_file>>tmp_str;
        }
        molpro_file>>tmp_str;
        molpro_file>>tmp_str;
        
        if (tmp_str=="active")
        {
            molpro_file>>tmp_str;
            molpro_file>>integer;
            *n_occ=integer;
            //cout<<" Number of active orbitals = "<<n_occ<<endl;//DEBOGAGE
            test=1;
        }
    }while(test!=1);//Search for the number of active orbitals in the molpro output file
    
    molpro_file.close();

    return 0;
}

int overlap_MO(double matrix[],int* n_occ,int* basis_size,std::string molpro_out_path)
{
    bool test(0);
    double *overlap;
    double *matrix2;
    double *MO_coeff_neutral;
    double *MO_coeff_cation;
    double floating;
    double *temp;
    double *temp2;
    
    using namespace std;
    
    string tmp_str;
    
    
    overlap=new double[*basis_size**basis_size];
    matrix2=new double [*n_occ**n_occ];
    
    
    
    ifstream molpro_file(molpro_out_path.c_str());
    
    if (!molpro_file.is_open())
    {
        cout<<"ERROR DURING MOLPRO OUTPUT FILE ACQUISITION"<<endl<<molpro_out_path.c_str()<<endl;
        return -1;
    }
    
    
    test=0;
    
    do
    {
        while (tmp_str!="MATRIX")
        {
            molpro_file>>tmp_str;
        }
        molpro_file>>tmp_str;
        
        if (tmp_str=="S")
        {
            molpro_file>>tmp_str;
            molpro_file>>tmp_str;
            molpro_file>>tmp_str;
            
            for (int i=0; i!=*basis_size**basis_size; i++)
            {
                
                 molpro_file>>floating;
                 overlap[i]=floating;
                //cout<<overlap[i]<<endl;//DEBOGAGE
            }
            test=1;
        }
    }while(test!=1);//Search for the overlap matrix in molpro output and save it to the overlap array
    
    test=0;
    
    molpro_file.seekg(ios_base::beg);
    
       MO_coeff_neutral=new double[*n_occ**basis_size];
    MO_coeff_cation=new double[*n_occ**basis_size];
    
    test=0;
    
    
    do
    {
        while (tmp_str!="NATURAL")
        {
            molpro_file>>tmp_str;
            if(molpro_file.eof())
                return -1;
        }
        
            molpro_file>>tmp_str;
        
        
        
        if (tmp_str=="ORBITALS")
        {
            //cout<<"Coefficients of the neutral"<<endl;//DEBOGAGE
            for (int i=0; i!=7+2**basis_size; i++)
            {
                molpro_file>>tmp_str;
                //cout<<tmp_str<<endl;//DEBOGAGE
            }
            
            for (int j=0; j!=*n_occ; j++)
            {
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                for (int k=0; k!=*basis_size; k++)
                {
                    molpro_file>>floating;
                    MO_coeff_neutral[j**basis_size+k]=floating;
                    //cout<<MO_coeff_neutral[j**basis_size+k]<<"     ";//DEBOGAGE
                }//std::cout<<endl;//DEBOGAGE
            }//cout<<endl;

            test=1;
        }
    }while(test!=1);//Search for the coefficients of the AO in the MO's of the neutral
    
    test=0;
    
    do
    {
        while (tmp_str!="NATURAL")
        {
            molpro_file>>tmp_str;
            if(molpro_file.eof())
                return -1;
        }
            molpro_file>>tmp_str;
        
        if (tmp_str=="ORBITALS")
        {
            //cout<<"Coefficients of the cation"<<endl;//DEBOGAGE
            for (int i=0; i!=7+2**basis_size; i++)
            {
                molpro_file>>tmp_str;
            }
            for (int j=0; j!=*n_occ; j++)
            {
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                molpro_file>>tmp_str;
                for (int k=0; k!=*basis_size; k++)
                {
                    
                    molpro_file>>floating;
                    MO_coeff_cation[j**basis_size+k]=floating;
                    //cout<<MO_coeff_cation[j**basis_size+k]<<"    ";//DEBOGAGE
                }//cout<<endl;
            }

            test=1;
        }
    }while(test!=1);//Search for the coefficients of the AO in the MO's of the cation
    
        molpro_file.close();
   

    
    
     /*for (int j=0; j!=*n_occ; j++)
     {

     for (int k=0; k!=*basis_size; k++)
     {
         cout<<MO_coeff_cation[j**basis_size+k]<<"     ";//DEBOGAGE
     }std::cout<<endl;//DEBOGAGE
     }cout<<endl;
     *///DEBOGAGE

    temp=new double[*basis_size**n_occ];
    temp2=new double[*basis_size**n_occ];
    
    transpose(MO_coeff_neutral, temp2, *n_occ, *basis_size);
    

    matrix_product(temp, overlap, temp2, *basis_size, *basis_size, *n_occ);
    /*
    for (int i=0;i!=*basis_size; i++)
    {
        for (int j=0; j!=*n_occ; j++)
        {
            cout<<temp[i**n_occ+j]<<"    ";
        }cout<<endl;
    }cout<<endl;*///DEBOGAGE
    
    matrix_product(matrix2, MO_coeff_cation, temp, *n_occ, *basis_size, *n_occ);
    
    transpose(matrix2, matrix, *n_occ, *n_occ);


    
    
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

    delete overlap;
    
    return 0;
}

int n_states_reader(int *n_states_neut,int *n_states_cat,std::string file_address)
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
    *n_states_neut=tmp_str.at(tmp_str.length()-2)-'0';

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
    *n_states_cat=tmp_str.at(tmp_str.length()-2)-'0';
    
    std::cout<<" cat states " <<*n_states_cat<<std::endl;
    molpro_output.close();
    return 0;
}

int num_of_ci_reader(int n_states_neut,int n_states_cat,int *n_ci_neut,int *n_ci_cat,std::string file_address)
{
    int position(0);
    std::string tmp_str("");
    int counter(0);
    double floating;
    
    std::ifstream molpro_output;
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
    
    *n_ci_neut=(counter-1)/(n_states_neut+1);
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
    
    *n_ci_cat=(counter-1)/(n_states_cat+1);
    
    
    molpro_output.close();
    
    return 0;
}

int ci_vec_reader(int n_states_neut,int n_states_cat,int n_occ,int n_elec_neut,int ci_size_neut,int ci_size_cat,double **ci_vector_neut,double **ci_vector_cat,std::string file_address)
{
    int position(0);
    std::string tmp_str("");
    int elec_index(0);
    double floating;
    
    std::ifstream molpro_output;
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;

    for (int i=0; i!=ci_size_neut; i++)
    {
        molpro_output>>tmp_str;
        for (int j=0; j!=n_occ; j++)
        {
            if(tmp_str.at(j)=='0')
                continue;
            
            else if(tmp_str.at(j)=='2')
            {
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=0;
                elec_index++;
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=1;
            }
            
            else if(tmp_str.at(j)=='a')
            {
                ci_vector_neut[0][(n_elec_neut+n_states_neut)*i+elec_index]=j;
                ci_vector_neut[1][(n_elec_neut)*i+elec_index]=0;
            }
            
            else if(tmp_str.at(j)=='b')
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
    molpro_output.close();
    
    
    if(!search(&position, file_address,position, "CI", 0, "vector"))
    {
        std::cout<<"CI VECTORS NOT FOUND IN MOLPRO OUPTUT FILE"<<std::endl;
        return 1;
    }
    
    molpro_output.open(file_address.c_str());
    
    molpro_output.seekg(position);
    
    molpro_output>>tmp_str;
    
    for (int i=0; i!=ci_size_cat; i++)
    {
        molpro_output>>tmp_str;
        for (int j=0; j!=n_occ; j++)
        {
            if(tmp_str.at(j)=='0')
                continue;
            
            else if(tmp_str.at(j)=='2')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=0;
                elec_index++;
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=1;
            }
            
            
            else if(tmp_str.at(j)=='a')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=0;
            }
            
            else if(tmp_str.at(j)=='b')
            {
                ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+elec_index]=j;
                ci_vector_cat[1][(n_elec_neut-1)*i+elec_index]=1;
            }
            
            elec_index++;
        }
        
        /*for (int v=0;v!=n_elec_neut-1 ; v++)
        {
            std::cout<<ci_vector_cat[0][(n_elec_neut-1+n_states_cat)*i+v]<<std::endl;
        }*/
        elec_index=0;
        
        for (int j=0; j!=n_states_cat; j++)
        {
            molpro_output>>floating;
            ci_vector_cat[0][(n_elec_neut-1+n_states_neut)*i+n_elec_neut-1+j]=floating;
            //std::cout<<std::stod(tmp_str)<<"   from config "<<i<<"   state "<<j<<std::endl;
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
    return 1;
}


