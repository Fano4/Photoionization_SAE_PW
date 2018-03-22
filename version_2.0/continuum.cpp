#include "continuum.hpp"


bool eikr(double *val,double k,double theta,double phi,double x,double y,double z)
{
   //This function returns the value of the plane wave defined by the spherical coordinates k,theta,phi at the point x,y,z using the vector val. val represents a complex number whose real part is val[0] and imaginary part is val[1]

   val[0]=cos((k*sin(theta)*cos(phi)*x + k*sin(theta)*sin(phi)*y+k*cos(theta)*z));//real part of the plane wave
   val[1]=sin((k*sin(theta)*cos(phi)*x + k*sin(theta)*sin(phi)*y+k*cos(theta)*z));//Imaginary part of the plane wave

   return 0;
}

bool pw_builder(double *Reinput,double *Iminput,double k,double *theta,double *phi,int angle_vec_size,double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz,int n_occ,double **neut_mo_cube_array)
{

   //This function builds a plane wave inside the pointers Reinput and Iminput
   double x;
   double y;
   double z;
   double Pi(acos(-1));
   double pw_value[2];

   using namespace std;


#pragma omp parallel for
   for( int i=0;i<angle_vec_size;i++)
   {
            for( int m=0;m<nx;m++)
            {
               x=xmin+m*(xmax-xmin)/nx;
               for( int n=0;n<ny;n++)
               {
                  y=ymin+n*(ymax-ymin)/ny;
                  for( int o=0;o<nz;o++)
                  {
                     z=zmin+o*(zmax-zmin)/nz;
                     eikr(pw_value, k,theta[i], phi[i], x, y, z);
                     Reinput[i*nx*ny*nz+m*ny*nz+n*nz+o]=pw_value[0]/sqrt((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
                     Iminput[i*nx*ny*nz+m*ny*nz+n*nz+o]=pw_value[1]/sqrt((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
                  }
               }
            }
   }
  // double *Retemp=new double[angle_vec_size];
  // double *Imtemp=new double[angle_vec_size];

/*   for(int i=0;i!=angle_vec_size;i++)
   {
      cube_dot_product(Reinput,&Reinput[i*nx*ny*nz],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,angle_vec_size,Retemp);
      cube_dot_product(Iminput,&Iminput[i*nx*ny*nz],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,angle_vec_size,Imtemp);
      std::cout<<"vector "<<i<<" =>"<<Retemp[i]+Imtemp[i]<<std::endl; 
   }*/
            pw_orthogonalizer(Reinput,Iminput,angle_vec_size,nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,n_occ,neut_mo_cube_array);


            return 1;
}
//double kmin,double kmax,int nk,int ntheta,int nphi

bool pw_orthogonalizer(double *Reinput,double *Iminput,int angle_vec_size,int nx,int ny,int nz,double dx,double dy,double dz,int n_occ,double **neut_mo_cube_array)
{

   const MKL_INT  num(nx*ny*nz);
   const int inc(1);
   const double coeff(1);
   double Redot[angle_vec_size];
   double Imdot[angle_vec_size];

//   std::cout<<"orthogonalizing..."<<std::endl;

   for(int k=0;k!=n_occ;k++)
   {
      cube_dot_product(Reinput,neut_mo_cube_array[k],nx,ny,nz,dx,dy,dz,angle_vec_size,Redot);
      cube_dot_product(Iminput,neut_mo_cube_array[k],nx,ny,nz,dx,dy,dz,angle_vec_size,Imdot);
/*      for( int i=0;i!=angle_vec_size;i++)
      {
         std::cout<<"dot prod value"<<Redot[i]<<" + i * "<<Imdot[i]<<std::endl;
      }*/

//      std::cout<<"updating basis..."<<std::endl;


#pragma omp parallel for
      for(int i=0;i<angle_vec_size;i++)
      {
         for( int m=0;m!=nx*ny*nz;m++)
         {
               Reinput[i*nx*ny*nz+m]-=Redot[i]*neut_mo_cube_array[k][m];
               Iminput[i*nx*ny*nz+m]-=Imdot[i]*neut_mo_cube_array[k][m];
         }
        // Redot[i]=-Redot[i];
        // Imdot[i]=-Imdot[i];
        // daxpy(&num,&Redot[i],neut_mo_cube_array[k],&inc,temp,&inc);
        // daxpy(&num,&coeff,temp,&inc,&Reinput[i*num],&inc);
        // dcopy(&num,temp2,&inc,temp,&inc);
        // daxpy(&num,&Imdot[i],neut_mo_cube_array[k],&inc,temp,&inc);
        // daxpy(&num,&coeff,temp,&inc,&Iminput[i*num],&inc);
        // dcopy(&num,temp2,&inc,temp,&inc);
      }

/*      for(int i=0;i!=nx;i++)
      {
         for(int j=0;j!=ny;j++)
         {
            for(int l=0;l!=nz;l++)
            {
               Reinput[i*ny*nz+j*nz+l]-=Redot*neut_mo_cube_array[k][i*ny*nz+j*nz+l];
               Iminput[i*ny*nz+j*nz+l]-=Imdot*neut_mo_cube_array[k][i*ny*nz+j*nz+l];
            }
         }
      }
*/
   }

//   std::cout<<"Done! ... ";
   return 1;
}

bool projector_lz(double *Reinput,double *Iminput,int angle_vec_size,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,int lz)
{
   double value(0);
   double phi(0);
   double x;
   double y;
   double dx=((xmax-xmin)/nx);
   double dy=((ymax-ymin)/ny);
   double *Reoutput=new double[nx*ny*nz];
   double *Imoutput=new double[nx*ny*nz];
   std::complex<double> output;

#pragma omp parallel for
   for(int n=0;n<angle_vec_size;n++)
   {
     for(int k=0;k!=nz;k++)
     {
        for(int i=0;i!=nx;i++)
        {
           x=xmin+i*dx;
           output=std::complex<double>(0,0);
           for(int j=0;j!=ny;j++)
           {
              y=ymin+j*dy;
              phi=atan2(x,y);
              output+=std::complex<double>(cos(lz*phi),-sin(lz*phi))*std::complex<double>(Reinput[n*nx*ny*nz+i*ny*nz+j*nz+k],Iminput[n*nx*ny*nz+i*ny*nz+j*nz+k])*(x/(pow(x,2)+pow(y,2)))*dy;
           }
           for(int j=0;j!=ny;j++)
           {
              y=ymin+j*dy;
              phi=atan2(x,y);
              Reoutput[i*ny*nz+j*nz+k]+=real(output*std::complex<double>(cos(lz*phi),sin(lz*phi)));
              Imoutput[i*ny*nz+j*nz+k]+=imag(output*std::complex<double>(cos(lz*phi),sin(lz*phi)));
           }
        }
        for(int j=0;j!=ny;j++)
        {
           y=ymin+j*dy;
           output=std::complex<double>(0,0);
           for(int i=0;i!=nx;i++)
           {
               x=xmin+i*dx;
               phi=atan2(x,y);
               output-=std::complex<double>(cos(lz*phi),-sin(lz*phi))*std::complex<double>(Reinput[n*nx*ny*nz+i*ny*nz+j*nz+k],Iminput[n*nx*ny*nz+i*ny*nz+j*nz+k])*(y/(pow(x,2)+pow(y,2)))*dx;
           }
           for(int i=0;i!=nx;i++)
           {
              x=xmin+i*dx;
              phi=atan2(x,y);
              Reoutput[i*ny*nz+j*nz+k]+=real(output*std::complex<double>(cos(lz*phi),sin(lz*phi)));
              Imoutput[i*ny*nz+j*nz+k]+=imag(output*std::complex<double>(cos(lz*phi),sin(lz*phi)));
              Reinput[n*nx*ny*nz+i*ny*nz+j*nz+k]=Reoutput[i*ny*nz+j*nz+k];
              Iminput[n*nx*ny*nz+i*ny*nz+j*nz+k]=Imoutput[i*ny*nz+j*nz+k];
           }
        }
     }
   }

   delete [] Reoutput;
   delete [] Imoutput;

   return 1;
}
