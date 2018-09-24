//
//  dyson_cube_writer.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 18/01/17.
//  Copyright © 2017 Stephan van den Wildenberg. All rights reserved.
//

#include "dyson_cube_writer.hpp"

bool cube_header(double *dyson_MO_basis_coeff,int n_occ,int n_states_neut,int n_states_cat,double **MO_cube_array,std::string Dyson_cube_loc,int state_neut,int state_cat,int n_nucl,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double *dyson_cube)
{
    using namespace std;
    double temp_norm(0);
    
    stringstream filename_str;
    string filename;
   /* for(int i=0;i!=n_states_neut*n_states_cat;i++)
    {
       temp_norm=0;  
       for(int v=0;v!=nx;v++)
       {
          for(int w=0;w!=ny;w++)
          {
             for(int y=0;y!=nz;y++)
             {
                temp_norm+=MO_cube_array[i][v*ny*nz+w*nz+y]*MO_cube_array[i][v*ny*ny+w*nz+y]*((xmax-xmin)/nx)*((ymax-ymin)/ny)*((zmax-zmin)/nz);
             }
          }
       }
       std::cout<<"dyson norm is "<<temp_norm<<std::endl;
       
       for(int v=0;v!=nx;v++)
       {
          for(int w=0;w!=ny;w++)
          {
             for(int y=0;y!=nz;y++)
             {
                MO_cube_array[i][v*ny*nz+w*nz+y]/=sqrt(temp_norm);
             }
          }
       }
    }
*/
   /* for(int k=0;k!=n_occ;k++)
    {
       std::cout<<dyson_MO_basis_coeff[n_states_cat*n_occ*state_neut+n_occ*state_cat+k]<<std::endl;
    }*/
    filename_str<<Dyson_cube_loc<<state_neut<<"."<<state_cat<<".cube";
    filename=filename_str.str();
    std::cout<<filename<<std::endl;
    
    ofstream Dyson_cube(filename.c_str());
    string temp;
    bool test(0);
    double floating;
    double sum;
    double temp_cube[nx*ny*nz];
    int index(0);
    int index2(0);
    int mo_index(0);

    Dyson_cube<<"Dyson orbital cube file"<<std::endl<<"Dyson between state "<<state_neut<<" of the neutral and state "<<state_cat<<" of the cation"<<"\n";
    Dyson_cube<<"-"<<n_nucl<<"  "<<xmin<<"  "<<ymin<<"  "<<zmin<<" \n";
    Dyson_cube<<nx<<"   "<<(xmax-xmin)/nx<<"    0.000000    0.000000 \n";
    Dyson_cube<<ny<<"   "<<"0.000000    "<<(ymax-ymin)/ny<<"    0.000000 \n";
    Dyson_cube<<nz<<"   "<<"0.000000    0.000000    "<<(zmax-zmin)/nz<<" \n";
    Dyson_cube<<"    1    1.000000    0.0000000000        0.0000000000      -1.295266 \n";
    Dyson_cube<<"    3    3.000000    0.0000000000        0.0000000000      0.217021 \n";
    //Dyson_cube<<"    6    6.000000    0.000000    0.993566    0.000000 \n";
    //Dyson_cube<<"    1    1.000000    1.775571    2.096624    0.000000 \n";
    //Dyson_cube<<"    1    1.000000   -1.775571    2.096624    0.000000 \n";
    //Dyson_cube<<"    8    8.000000    0.000000   -1.269371    0.000000 \n";
    //Dyson_cube<<"  6    6.000000    0.000000    0.000000    1.173348 \n";
    //Dyson_cube<<"  7    7.000000    0.000000    0.000000   -1.005506 \n";
    //Dyson_cube<<"1    1.000000    0.000000    0.000000    3.212534 \n";

    Dyson_cube<<"1  111 \n";

        sum=0;
           for(int i=0;i!=nx;i++)
           {
              for(int j=0;j!=ny;j++)
              {
                 for(int l=0;l!=nz;l++)
                 {
                    dyson_cube[i*ny*nz+j*nz+l]=0;
                 }
              }
           }
        for (int k=0; k!=n_occ; k++)
        {
           for(int i=0;i!=nx;i++)
           {
              for(int j=0;j!=ny;j++)
              {
                 for(int l=0;l!=nz;l++)
                 {
                    dyson_cube[i*ny*nz+j*nz+l]+=dyson_MO_basis_coeff[n_occ*n_states_cat*state_neut+n_occ*state_cat+k]*MO_cube_array[k][i*ny*nz+j*nz+l];
                 }
              }
           }
        }
        temp_norm=0;
        for(int i=0;i!=nx;i++)
        {
           for(int j=0;j!=ny;j++)
           {
              for(int k=0;k!=nz;k++)
              {
                Dyson_cube<<scientific<<setw(16)<<dyson_cube[i*ny*nz+j*nz+k];
                temp_norm+=dyson_cube[i*ny*nz+j*nz+k]*dyson_cube[i*ny*nz+j*nz+k]*((xmax-xmin)/nx)*((ymax-ymin)/ny)*((zmax-zmin)/nz);
                index++;
                index2++;
                if (index2%6==0)
                {
                    Dyson_cube<<"\n";
                }
                if (index%ny==0)
                {
                    index2=0;
                    Dyson_cube<<"\n";
                }
              }
           }
        }
        std::cout<<"Dyson_norm is "<<temp_norm<<std::endl;
        
        Dyson_cube.close();
    
    return 1;
}

bool cube_reader(int mo_index1,int mo_index2,int nx,int ny,int nz,std::string MO_cube_loc,double *cube_array)
{
   using namespace std;
   string temp;
   int num_of_nucl(0);
   string filename;
   stringstream filename_str;
   stringstream strnum_str;

        filename_str.str("");
        filename_str<<MO_cube_loc.c_str()<<mo_index1+1<<"."<<mo_index2+1<<".cube";
        filename=filename_str.str();

        ifstream MO_cube_out;
        
        MO_cube_out.open(filename.c_str());
         if (!MO_cube_out.is_open())
         {
             cout<<"ERROR: CANNOT OPEN CUBE FILE "<<filename.c_str()<<endl;
             return 0;
         }
         MO_cube_out>>temp;
         MO_cube_out>>temp;
         MO_cube_out>>temp;
         num_of_nucl=fabs(atoi(temp.c_str()));
         //std::cout<<"read number of atoms: "<<num_of_nucl<<std::endl;

         for(int i=0;i!=15+fabs(num_of_nucl)*5+2;i++)
         {
            MO_cube_out>>temp;
//            std::cout<<temp<<std::endl;
         }

        for(int i=0;i!=nx*ny*nz;i++)
        {
           cube_array[i]=0;
            MO_cube_out>>cube_array[i];
        }
        MO_cube_out.close();

   return 1;
}
