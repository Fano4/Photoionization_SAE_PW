//
//  dyson_cube_writer.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 18/01/17.
//  Copyright © 2017 Stephan van den Wildenberg. All rights reserved.
//

#include "dyson_cube_writer.hpp"

bool cube_header(double *dyson_MO_basis_coeff,int n_occ,int n_states_neut,int n_states_cat,std::string MO_cube_loc,std::string Dyson_cube_loc)
{
    
    using namespace std;
    
    stringstream filename_str;
    string filename;
    ifstream MO_cube_out;
    ofstream Dyson_cube;
    ifstream Dyson_cube_read;
    string temp;
    bool test(0);
    char cursor;
    int match_loc(0);
    int temp_loc(0);
    double sum_vector[n_occ];
    double floating;
    double sum;
    int i(0);
    int j(0);
    int num_of_points_cube(80);
    int index(0);
    int index2(0);
    
    std::cout<<"CHECK 1"<<std::endl;
    
    for (int k=0; k!=n_occ; k++)
    {
        filename_str.str("");
        filename_str<<MO_cube_loc.c_str()<<k+1<<".1.cube";
    std::cout<<"CHECK 2"<<std::endl;
        filename=filename_str.str();
        
        MO_cube_out.open(filename.c_str());
         if (!MO_cube_out.is_open())
         {
             cout<<"ERROR: CANNOT OPEN MOLPRO CUBE FILE "<<filename.c_str()<<endl;
             return 0;
         }
        else
        {
            //cout<<filename.c_str()<<" file open "<<endl;
            MO_cube_out.close();
        }
    }
    
    filename_str.str("");
    filename_str<<MO_cube_loc.c_str()<<"1.1.cube";
    filename=filename_str.str();
    
    
    search(&match_loc, filename, 0, "1",0,"101");
    
    MO_cube_out.open(filename.c_str());
    
    filename_str.str("");
    filename_str<<Dyson_cube_loc.c_str()<<i<<"."<<j<<".cube";
    filename=filename_str.str();
    
    Dyson_cube.open(filename.c_str());
    Dyson_cube.seekp(ios_base::beg);
    
    do
    {
        cursor=MO_cube_out.get();
        Dyson_cube<<cursor;
    }while(int(MO_cube_out.tellg())!=match_loc);
    
    Dyson_cube<<"\n";
    Dyson_cube.close();
    MO_cube_out.close();
    
    filename_str.str("");
    filename_str<<Dyson_cube_loc.c_str()<<i<<"."<<j<<".cube";
    filename=filename_str.str();
    
    
    do
    {
        for (int k=0; k!=n_occ; k++)
        {
            filename_str.str("");
            filename_str<<MO_cube_loc.c_str()<<k+1<<".1.cube";
            filename=filename_str.str();
        
            MO_cube_out.open(filename.c_str());
            MO_cube_out.seekg(match_loc);
            MO_cube_out>>floating;
            sum_vector[k]=floating;
            
            test=MO_cube_out.eof();
            temp_loc=MO_cube_out.tellg();
            
            MO_cube_out.close();
        }
        sum=0;
        for (int k=0; k!=n_occ; k++)
        {
            sum+=dyson_MO_basis_coeff[n_occ*n_states_cat*i+n_occ*j+k]*sum_vector[k];
        }
        std::cout<<index<<endl;//DEBOGAGE
        filename_str.str("");
        filename_str<<Dyson_cube_loc.c_str()<<i<<"."<<j<<".cube";
        filename=filename_str.str();
        Dyson_cube.open(filename.c_str());
        Dyson_cube.seekp(match_loc);
        match_loc=temp_loc;
        Dyson_cube<<scientific<<setw(16)<<sum;
        index++;
        index2++;
        if (index2%6==0)
        {
            Dyson_cube<<"\n";
        }
        if (index%num_of_points_cube==0)
        {
            index2=0;
            Dyson_cube<<"\n";
        }
        
        Dyson_cube.close();
        
    }while (!test);

    
    return 1;
}

bool cube_reader(int mo_index,int nx,int ny,int nz,std::string MO_cube_loc,double *cube_array)
{
   using namespace std;
   ifstream MO_cube_out;
   string strnum;
   string filename;
   stringstream filename_str;
   stringstream strnum_str;
   int match_loc;

        filename_str.str("");
        filename_str<<MO_cube_loc.c_str()<<mo_index+1<<".1.cube";
        filename=filename_str.str();
        
        MO_cube_out.open(filename.c_str());
         if (!MO_cube_out.is_open())
         {
             cout<<"ERROR: CANNOT OPEN MOLPRO CUBE FILE "<<filename.c_str()<<endl;
             return 0;
         }
        else
        {
            //cout<<filename.c_str()<<" file open "<<endl;
            MO_cube_out.close();
        }
        strnum_str.str("");
        strnum_str<<(mo_index+1)<<"01";
        strnum=strnum_str.str();

        search(&match_loc, filename, 0, "1",0,strnum.c_str());

        MO_cube_out.open(filename.c_str());

        MO_cube_out.seekg(match_loc);

        for(int i=0;i!=nx;i++)
        {
           for(int j=0;j!=ny;j++)
           {
              for(int k=0;k!=nz;k++)
              {
                 MO_cube_out>>cube_array[i*ny*nz+j*nz+k];
              }
           }
        }

   return 1;
}
