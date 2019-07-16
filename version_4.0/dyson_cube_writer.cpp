//
//  dyson_cube_writer.cpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 18/01/17.
//  Copyright © 2017 Stephan van den Wildenberg. All rights reserved.
//

#include "dyson_cube_writer.hpp"

bool cube_header(int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int *n_nucl_dim,int *grid_size,int *num_of_nucl,int* basis_size,int *contraction_number,double nucl_coord,double **nucl_spher_pos,double *MO_coeff_neutral,double *dyson_mo_coeff,double **contraction_coeff,double **contraction_zeta,int* nucleus_basis_func,std::string *basis_func_type,int **angular_mom_numbers,std::string Dyson_cube_loc,int state_neut,int state_cat,int n_nucl,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax)
{
    using namespace std;
    double temp_norm(0);

    double* dyson_cube=new double[nx*ny*nz];
   double r(0);
   double thet(0);
   double phi(0);
   double x(0);
   double y(0);
   double z(0);
    
    stringstream filename_str;
    string filename;

    filename_str<<Dyson_cube_loc<<state_neut<<"."<<state_cat<<".cube";
    filename=filename_str.str();
    std::cout<<filename<<std::endl;
    
    ofstream Dyson_cube;
    string temp;
    bool test(0);
    double floating;
    double temp_cube[nx*ny*nz];
    int index(0);
    int index2(0);
    int mo_index(0);
    int i(0);
    int j(0);
    int k(0);
    int l(0);

           for( i=0;i!=nx;i++)
           {
              for( j=0;j!=ny;j++)
              {
                 for( l=0;l!=nz;l++)
                 {
                    dyson_cube[i*ny*nz+j*nz+l]=0;
                 }
              }
           }
      for( i=0;i<nx;i++)
      {
         std::cout<<"Writing cube "<<i<<std::endl;
         x=xmin+i*(xmax-xmin)/nx;
#pragma omp parallel for private(j,k,l,r,thet,phi,y,z) shared(i,x,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,dyson_cube,dyson_mo_coeff,nucl_spher_pos,nucleus_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size,n_occ,angular_mom_numbers)
         for( j=0;j<ny;j++)
         {
            y=ymin+j*(ymax-ymin)/ny;
            for( l=0;l<nz;l++)
            {
               z=zmin+l*(zmax-zmin)/nz;
                for (k=0; k<*n_occ; k++)
                {
                    dyson_cube[i*ny*nz+j*nz+l]+=dyson_mo_coeff[*n_states_cat**n_occ*state_neut+state_cat**n_occ+k]*MO_value(k,x,y,z,nucl_spher_pos,nucleus_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,*basis_size,angular_mom_numbers);
 //                   std::cout<<dyson_mo_coeff[*n_states_cat**n_occ*state_neut+state_cat**n_occ+k]<<"*"<<MO_value(k,x,y,z,nucl_spher_pos,nucleus_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,*basis_size,angular_mom_numbers)<<std::endl;
                 }
              }
           }
        }

      std::cout<<"cube built!"<<std::endl;
       temp_norm=0;  
       for( i=0;i!=nx;i++)
       {
          for( j=0;j!=ny;j++)
          {
             for( l=0;l!=nz;l++)
             {
                temp_norm+=dyson_cube[i*ny*nz+j*nz+l]*dyson_cube[i*ny*nz+j*nz+l]*((xmax-xmin)/nx)*((ymax-ymin)/ny)*((zmax-zmin)/nz);
             }
          }
       }
       std::cout<<"dyson norm is "<<temp_norm<<std::endl;

    Dyson_cube.open(filename.c_str());
    Dyson_cube<<"Dyson orbital cube file"<<std::endl<<"Dyson between state "<<state_neut<<" of the neutral and state "<<state_cat<<" of the cation"<<"\n";
    Dyson_cube<<"-"<<n_nucl<<"  "<<xmin<<"  "<<ymin<<"  "<<zmin<<" \n";
    Dyson_cube<<nx<<"   "<<(xmax-xmin)/nx<<"    0.000000    0.000000 \n";
    Dyson_cube<<ny<<"   "<<"0.000000    "<<(ymax-ymin)/ny<<"    0.000000 \n";
    Dyson_cube<<nz<<"   "<<"0.000000    0.000000    "<<(zmax-zmin)/nz<<" \n";
    Dyson_cube<<"    1    1.000000    0.0000000000        0.0000000000      "<<-1.295266/0.529<<" \n";
    Dyson_cube<<"    3    3.000000    0.0000000000        0.0000000000      "<<0.217021/0.529<<" \n";
    Dyson_cube<<"1  111 \n";
        for( i=0;i!=nx;i++)
        {
           for( j=0;j!=ny;j++)
           {
              for( l=0;l!=nz;l++)
              {
                Dyson_cube<<scientific<<setw(16)<<dyson_cube[i*ny*nz+j*nz+l];
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
            //std::cout<<temp<<std::endl;
         }

        for(int i=0;i!=nx*ny*nz;i++)
        {
            MO_cube_out>>cube_array[i];
        }
        MO_cube_out.close();

   return 1;
}

