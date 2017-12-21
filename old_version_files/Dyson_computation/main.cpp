//
//  main.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 19/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#include "dyson_computation.hpp"

int main(int argc, const char * argv[])
{
    
    std::string molpro_out_path="/Users/Fano/Documents/Archivage ULG /Doctorat/Recherche/H2O/H2O_test_info.out";
    std::string MO_cube_path="/Users/Fano/Documents/Archivage ULG /Doctorat/Recherche/H2O/H2O_test/h2o_neut_orbital_";
    std::string Dyson_cube_path="/Users/Fano/Documents/Archivage ULG /Doctorat/Recherche/H2O/H2O_Dyson_";
    int n_occ;
    int basis_size;
    double *array;
    int n_states_neut(0);
    int n_states_cat(0);
    int n_elec_neut(10);
    int ci_size_neut(0);
    int ci_size_cat(0);
    double *ci_vec_neut[2];
    double *ci_vec_cat[2];
    double *temp;
    bool test(0);
    bool test2(0);
    
    double *Dyson_MO_basis_coeff;
    
    //int dim1=2;
    //int dim2=3;
    //int dim3=4;
    //double testA[81]={0.99993,0,-0.0086651,-0,0,0.0078478,0,-4.3368e-19,-0,0,0.99993,-0,-0.0086651,0,0,0.0078478,-0,-4.3368e-19,0.010464,0,0.96279,0,0,-0.2694,-0,0,0,-0.005215,-0,0.26933,0,0,0.96173,0,-3.4694e-18,-0,-0,-0.005215,0,0.26933,0,0,0.96173,-0,-3.4694e-18,-6.5052e-19,-0,-1.2143e-17,-0,0,0,0,0.9988,0,-0,-6.5052e-19,-0,-1.2143e-17,0,0,0,0,0.9988,0,0,0,0,-0.032007,0,0,0,0,0,0,0,0,-0,0,0,0,0};
    //double testB[12]={1,2,1,2,1,2,1,2,3,4,3,4};
    //double testC[12]={};
    
    //transpose(testB, testC, 3,4);
    
    //matrix_product(testC, testA, testB, dim1, dim2, dim3);
   /*
    for (int i=0; i!=dim1; i++)
    {
        for (int j=0; j!=dim3; j++)
        {
            std::cout<<testC[i*dim3+j]<<"    ";
        }std::cout<<std::endl;
    }std::cout<<std::endl;*/
    
    //std::cout<<determinant(testA, 9)<<std::endl;

  
    size_query(&n_occ, &basis_size, molpro_out_path);
    
    //std::cout<<n_occ<<"   "<<basis_size<<std::endl;
    
    array=new double[n_occ*n_occ];
    
    overlap_MO(&array[0],&n_occ,&basis_size,molpro_out_path);
    

    for (int i=0; i!=n_occ; i++)
    {
        for (int j=0; j!=n_occ; j++)
        {
            std::cout<<std::setw(12)<<std::setprecision(5)<<array[i*n_occ+j]<<"   ";
            
        }std::cout<<std::endl;
    }std::cout<<std::endl<<std::endl;
    
    n_states_reader(&n_states_neut, &n_states_cat, molpro_out_path);
    
    num_of_ci_reader(n_states_neut, n_states_cat, &ci_size_neut, &ci_size_cat, molpro_out_path);
    
    
    //std::cout<<ci_size_cat<<"    "<<ci_size_neut<<std::endl;
    
    ci_vec_neut[0]=new double[n_elec_neut*ci_size_neut+n_states_neut*ci_size_neut];//vector partitionned in two sections. Section 1 is filled with the mo label of each electron. section 2 is filled with CI coeff.
    ci_vec_neut[1]=new double[n_elec_neut*ci_size_neut];//this vector represents the spin state of each electron
    ci_vec_cat[0]=new double [(n_elec_neut-1)*ci_size_cat+n_states_cat*ci_size_cat];
    ci_vec_cat[1]=new double [(n_elec_neut-1)*ci_size_cat];
    Dyson_MO_basis_coeff=new double[n_occ*n_states_neut*n_states_cat];
    temp = new double[(n_elec_neut-1)*(n_elec_neut-1)];
    

    
    
    
    
    //std::cout<<pow(factorial(n_occ)/(factorial(n_elec_neut/2)*factorial(n_occ-n_elec_neut/2)),2)<<std::endl;
    
    ci_vec_reader(n_states_neut, n_states_cat, n_occ, n_elec_neut,ci_size_neut,ci_size_cat, ci_vec_neut, ci_vec_cat, molpro_out_path);
    
    /*for (int j=0; j!=ci_size_neut; j++)
    {
        for (int i=0; i!=n_occ; i++)
        {
            std::cout<<ci_vec_neut[0][(n_occ+n_states_neut)*j+i]<<"    ";
        }std::cout<<std::endl;
        for (int i=0; i!=n_states_neut; i++)
        {
            std::cout<<ci_vec_neut[0][(n_occ+n_states_neut)*j+n_occ+i]<<"    ";
        }std::cout<<std::endl;
        
    }*/
    
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
                        
                        for (int m=0; m!=n_elec_neut; m++)  //Over the electrons of the neutral
                        {
                            for (int p=0; p!=n_occ; p++)//Over the MO of the neutral
                            {
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
                                                    temp[(n_elec_neut-1)*(m-test2)+o]=array[n_occ*p+q]*kronecker_delta(ci_vec_neut[1][n_elec_neut*n+(m-test2)], ci_vec_cat[1][(n_elec_neut-1)*l+o]);
                                                }
                                            }
                                        }
                                        
                                    }
                                }

                            }
                            
                        }
                        
                        
                        /*//====================CHECKING DETERMINANT===============
                        
                        /for (int m=0; m!=n_elec_neut; m++)
                        {
                            std::cout<<(ci_vec_neut[1][n_elec_neut*n+m])<<std::endl;
                            
                        }std::cout<<std::endl;
                        for (int m=0; m!=n_elec_neut-1; m++)
                        {
                            std::cout<<(ci_vec_cat[1][(n_elec_neut-1)*l+m])<<std::endl;
                            
                        }std::cout<<std::endl;/
                        for (int m=0; m!=n_elec_neut-1; m++)
                        {
                            for (int o=0; o!=n_elec_neut-1; o++)
                            {
                                std::cout<<std::setw(12)<<std::setprecision(5)<<temp[(n_elec_neut-1)*m+o];
                            }std::cout<<std::endl;
                        }std::cout<<std::endl;
                        
                        
                        
                        *///====================CHECKING DETERMINANT===============
                        if(test2)
                        Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]+=ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]*ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]*determinant(temp,(n_elec_neut-1));

                        
                        //std::cout<<std::endl<<std::setw(12)<<std::setprecision(5)<<"config neutral "<<n<<" = "<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]<<std::setw(12)<<std::setprecision(5)<<"  config cation "<<l<<" = "<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]<<std::setw(12)<<std::setprecision(5)<<determinant(temp,(n_elec_neut-1))<<std::endl<<std::endl;
                        
                        //std::cout<<"config neutral "<<n<<" = "<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*n+n_elec_neut+i]<<"  config cation "<<l<<" = "<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*l+n_elec_neut-1+j]<<std::endl;
                    }
                }
                if(test)
                {
                    std::cout<<std::endl<<"states "<<i<<"  and  "<<j<<"    "<<Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]<<"   MO  "<<k<<std::endl<<"====================================="<<std::endl<<std::endl;
                }
                else
                {
                    Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]=0;
                    std::cout<<std::endl<<"states "<<i<<"  and  "<<j<<"    "<<Dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]<<"   MO  "<<k<<std::endl<<"====================================="<<std::endl<<std::endl;
                }
                    
            }
        }
    }
    
    cube_header(Dyson_MO_basis_coeff, n_occ, n_states_neut, n_states_cat, MO_cube_path,Dyson_cube_path);
    
    
    return 0;
    
}
