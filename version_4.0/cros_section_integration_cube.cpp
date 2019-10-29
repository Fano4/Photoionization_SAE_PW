#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <complex>

bool write_mfpad(double k,int n_theta,int n_phi,int nx,double xmin,double xmax,int ny,double ymin,double ymax,int nz,double zmin,double zmax,double *cs_cube,double *mfpad,std::string mfpad_adress="");
double total_cs(int n_theta,int n_phi,double **sphere_dist,double *mfpad);
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address);
void spherical_extract_from_cube(double k,double** sphere_dist,double n_points_sphere,double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz,double* cube,double* sphere);
void sphere_dist_gen(int n_points_sphere,double **sphere_dist,bool randiso,int n_phi);
void cube_diff_cross_section(double x_comp,double y_comp,double z_comp,int cube_size,double* cs_cube,int neut_state,int cat_state,std::string pice_address,double position_nuc);

int main(int argc, char *argv [])
{
   //THIS CODE COMPUTES THE PHOTOIONIZATION CROSS SECTION FROM CUBE FILES IN RECIPROCAL SPACE. FIRST STEP IS COMPUTATION OF CROSS SECTION IN THE CUBE FROM PICE. SECOND STEP IS CONVERTING THE CROSS SECTION FROM CARTESIAN TO SPHERICAL COORDINATES. THEN, EITHER MPFAD IS COMPUTED AS DIFFERENTIAL CROSS SECTION OR PHOTOELECTRON SPECTRUM IS COMPUTED BY INTEGRATION OF ALL MFPAD AS A FUNCTION OF K
   //
   using namespace std;
   string pice_address("/data1/home/stephan/LiH_512_points_pice/LiH_");
   string spectrum_address("/data1/home/stephan/mfpad_test.txt");
   int n_theta(90);
   int n_phi(120);
   int n_points_sphere(128);
   int n_points(n_theta*n_phi);
   double **sphere_dist=new double *[2];
   sphere_dist[0]=new double[n_points];
   sphere_dist[1]=new double[n_points];
   int nx(64);//125);
   int ny(64);//125);
   int nz(64);//125);
   double xmin(-acos(-1)/2);
   double xmax(acos(-1)/2);
   double ymin(-acos(-1)/2);
   double ymax(acos(-1)/2);
   double zmin(-acos(-1)/2);
   double zmax(acos(-1)/2);

   double kmax(acos(-1)/2);
   double nk=500;

   double theta_elec(0);//(acos(-1)/2);
   double phi_elec(0);//(acos(-1)/2);
   double x_comp(sin(theta_elec)*cos(phi_elec));
   double y_comp(sin(theta_elec)*sin(phi_elec));
   double z_comp(cos(theta_elec));
   double *mfpad=new double[n_points];
   double *cs_cube=new double[nx*ny*nz];
   const double Pi(acos(-1));

   double k(0.66407689747345189);

   std::cout<<"Electric field components : (X-comp) , (Y-comp) , (Z-comp)"<<std::endl<<x_comp<<" , "<<y_comp<<", "<<z_comp<<std::endl;

   cube_diff_cross_section(x_comp,y_comp,z_comp,nx*ny*nz,cs_cube,3,0,pice_address,1.6125);
   //write_mfpad(k,n_theta,n_phi,nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax,cs_cube,mfpad,spectrum_address);
   sphere_dist_gen(n_points,sphere_dist,0,n_phi);
   for(int i=1;i!=nk;i++)
   {
      k=i*kmax/nk;
      write_mfpad(k,n_theta,n_phi,nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax,cs_cube,mfpad);
      std::cout<<k<<"   "<<k*k*total_cs(n_theta,n_phi,sphere_dist,mfpad)<<std::endl;;
   }

   return 0;
}

bool write_mfpad(double k,int n_theta,int n_phi,int nx,double xmin,double xmax,int ny,double ymin,double ymax,int nz,double zmin,double zmax,double *cs_cube,double *mfpad,std::string mfpad_adress)
{
   using namespace std;
   int n_points_sphere=n_theta*n_phi;
   double **sphere_dist=new double*[2];
   sphere_dist[0]=new double[n_points_sphere];
   sphere_dist[1]=new double[n_points_sphere];

//   std::cout<<"generating spherical distribution"<<std::endl;
   sphere_dist_gen(n_points_sphere,sphere_dist,0,n_phi);
//   std::cout<<"extracting sphere from cube"<<std::endl;
   spherical_extract_from_cube(k,sphere_dist,n_points_sphere,xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,cs_cube,mfpad);

//   std::cout<<"writing mfpad"<<std::endl;

   if(mfpad_adress!="")
   {
      ofstream mfpad_out;
      mfpad_out.open(mfpad_adress.c_str());
      for(int i=0;i!=n_theta;i++)
      {
         for(int j=0;j!=n_phi;j++)
         {
            mfpad_out<<k<<"   "<<sphere_dist[0][i*n_phi+j]<<"   "<<sphere_dist[1][i*n_phi+j]<<"   "<<mfpad[i*n_phi+j]<<std::endl;
         }mfpad_out<<std::endl;
      }
      mfpad_out.close();
   }
   return 1;
}
double total_cs(int n_theta,int n_phi,double **sphere_dist,double *mfpad)
{
   double sum(0);
   double Pi=acos(-1);
   for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
          sum+=sin(sphere_dist[0][i*n_phi+j])*(Pi/n_theta)*(2*Pi/n_phi)*mfpad[i*n_phi+j];
      }
   }
   return sum;
}
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address)
{
   double * k=new double[n_points];
   double * theta=new double[n_points];
   double * phi=new double[n_points];
   double *RePICE_x=new double[n_points];
   double *ImPICE_x=new double[n_points];
   double *RePICE_y=new double[n_points];
   double *ImPICE_y=new double[n_points];
   double *RePICE_z=new double[n_points];
   double *ImPICE_z=new double[n_points];

   const double pi(acos(-1));
   const double Lx(26.6337);
   const double Ly(26.6337);
   const double Lz(27.876);
   const int nk(1000);
   double *cs=new double[nk];

   using namespace std;
   ifstream input;
   input.open(dipole_address.c_str());

   if(!input.is_open())
   {      
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=n_points;i++)
      {
         input>>k[i];
         input>>theta[i];
         input>>phi[i];
         input>>RePICE_x[i];
         input>>ImPICE_x[i];
         input>>RePICE_y[i];
         input>>ImPICE_y[i];
         input>>RePICE_z[i];
         input>>ImPICE_z[i];
      }
      for(int i=0;i!=nk;i++)
      {
         cs[i]=0;
         for(int j=0;j!=n_points/nk;j++)
         {
            cs[i]+=(4*pi*nk/n_points)*(2*k[i*n_points/nk+j]*Lx*Ly*Lz/(137.0396096*pi))*real(std::complex<double>(x_comp*RePICE_x[i*n_points/nk+j]+y_comp*RePICE_y[i*n_points/nk+j]+z_comp*RePICE_z[i*n_points/nk+j],x_comp*ImPICE_x[i*n_points/nk+j]+y_comp*ImPICE_y[i*n_points/nk+j]+z_comp*ImPICE_z[i*n_points/nk+j])*std::complex<double>(x_comp*RePICE_x[i*n_points/nk+j]+y_comp*RePICE_y[i*n_points/nk+j]+z_comp*RePICE_z[i*n_points/nk+j],-(x_comp*ImPICE_x[i*n_points/nk+j]+y_comp*ImPICE_y[i*n_points/nk+j]+z_comp*ImPICE_z[i*n_points/nk+j])));
         }
         std::cout<<pow(k[i*n_points/nk],2)*27.211/2<<", "<<cs[i]<<std::endl;
      }
      input.close(); 
   }
   return 1;
}
void spherical_extract_from_cube(double k,double** sphere_dist,double n_points_sphere,double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz,double* cube,double* sphere)
{
   double xsphere(0);
   double ysphere(0);
   double zsphere(0);
   int x_index(0);
   int y_index(0);
   int z_index(0);

   for(int i=0;i!=n_points_sphere;i++)
   {
      xsphere=k*sin(sphere_dist[0][i])*cos(sphere_dist[1][i]);
      ysphere=k*sin(sphere_dist[0][i])*sin(sphere_dist[1][i]);
      zsphere=k*cos(sphere_dist[0][i]);

      x_index=int(round((xsphere-xmin)*nx/(xmax-xmin)));
      y_index=int(round((ysphere-ymin)*ny/(ymax-ymin)));
      z_index=int(round((zsphere-zmin)*nz/(zmax-zmin)));

      sphere[i]=k*cube[x_index*ny*nz+y_index*nz+z_index];
   }
}
void cube_diff_cross_section(double x_comp,double y_comp,double z_comp,int cube_size,double* cs_cube,int neut_state,int cat_state,std::string pice_address,double position_nuc)
{
   using namespace std;

   string test_str;

   ifstream repicex;
   ifstream impicex;
   ifstream repicey;
   ifstream impicey;
   ifstream repicez;
   ifstream impicez;
   string s_repice_x;
   stringstream ss_repice_x;
   string s_impice_x;
   stringstream ss_impice_x;
   string s_repice_y;
   stringstream ss_repice_y;
   string s_impice_y;
   stringstream ss_impice_y;
   string s_repice_z;
   stringstream ss_repice_z;
   string s_impice_z;
   stringstream ss_impice_z;
   bool test(0);
   ss_repice_x.str("");
   ss_repice_x<<pice_address.c_str()<<"RePICE_"<<position_nuc<<"_X_"<<neut_state<<"_"<<cat_state<<".txt";
   s_repice_x=ss_repice_x.str();
   ss_impice_x.str("");
   ss_impice_x<<pice_address.c_str()<<"ImPICE_"<<position_nuc<<"_X_"<<neut_state<<"_"<<cat_state<<".txt";
   s_impice_x=ss_impice_x.str();
   ss_repice_y.str("");
   ss_repice_y<<pice_address.c_str()<<"RePICE_"<<position_nuc<<"_Y_"<<neut_state<<"_"<<cat_state<<".txt";
   s_repice_y=ss_repice_y.str();
   ss_impice_y.str("");
   ss_impice_y<<pice_address.c_str()<<"ImPICE_"<<position_nuc<<"_Y_"<<neut_state<<"_"<<cat_state<<".txt";
   s_impice_y=ss_impice_y.str();
   ss_repice_z.str("");
   ss_repice_z<<pice_address.c_str()<<"RePICE_"<<position_nuc<<"_Z_"<<neut_state<<"_"<<cat_state<<".txt";
   s_repice_z=ss_repice_z.str();
   ss_impice_z.str("");
   ss_impice_z<<pice_address.c_str()<<"ImPICE_"<<position_nuc<<"_Z_"<<neut_state<<"_"<<cat_state<<".txt";
   s_impice_z=ss_impice_z.str();

   std::cout<<"Opening cube files"<<std::endl;

   repicex.open(s_repice_x.c_str());
   if(!repicex.is_open())
   {
      std::cout<<"Error while opening "<<s_repice_x.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      repicex>>test_str;
   }while(test_str != "111");
   impicex.open(s_impice_x.c_str());
   if(!impicex.is_open())
   {
      std::cout<<"Error while opening "<<s_impice_x.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      impicex>>test_str;
   }while(test_str != "111");
   repicey.open(s_repice_y.c_str());
   if(!repicey.is_open())
   {
      std::cout<<"Error while opening "<<s_repice_y.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      repicey>>test_str;
   }while(test_str != "111");
   impicey.open(s_impice_y.c_str());
   if(!impicey.is_open())
   {
      std::cout<<"Error while opening "<<s_impice_y.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      impicey>>test_str;
   }while(test_str != "111");
   repicez.open(s_repice_z.c_str());
   if(!repicez.is_open())
   {
      std::cout<<"Error while opening "<<s_repice_z.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      repicez>>test_str;
   }while(test_str != "111");
   impicez.open(s_impice_z.c_str());
   if(!impicez.is_open())
   {
      std::cout<<"Error while opening "<<s_impice_z.c_str()<<std::endl;
      exit(EXIT_FAILURE);
   }
   do
   {
      impicez>>test_str;
   }while(test_str != "111");

   test=(repicex.is_open())*(impicex.is_open())*(repicey.is_open())*(impicey.is_open())*(repicez.is_open())*(impicez.is_open());
   if(!test)
   {
      std::cout<<"FATAL ERROR WHILE OPENING PICE CUBE FILE"<<std::endl;
      exit(EXIT_FAILURE);
   }

   double repice_x;
   double impice_x;
   double repice_y;
   double impice_y;
   double repice_z;
   double impice_z;

   std::cout<<"extracting data from cube files"<<std::endl;
   for(int i=0;i!=cube_size;i++)
   {
      repicex>>repice_x;
      impicex>>impice_x;
      repicey>>repice_y;
      impicey>>impice_y;
      repicez>>repice_z;
      impicez>>impice_z;
      cs_cube[i]=(pow(128/(2*acos(-1)),3)*16*acos(-1)*acos(-1))*norm(complex<double>(x_comp*repice_x+y_comp*repice_y+z_comp*repice_z,x_comp*impice_x+y_comp*impice_y+z_comp*impice_z))/137;
   }
   std::cout<<"probe end of cube reading"<<std::endl;
   repicex.close();
   impicex.close();
   repicey.close();
   impicey.close();
   repicez.close();
   impicez.close();
}
void sphere_dist_gen(int n_points_sphere,double **sphere_dist,bool randiso=1,int n_phi=0)
{
   double random;
   double random2;
   double random3;
   double temp;
   double n_theta(0);
   double theta(0);
   double phi(0);
   double Pi=acos(-1);

   if(randiso)
   {
      for(int i=0;i!=n_points_sphere;i++)
      {
         random=double(rand()%2000)-1000;
         random2=double(rand()%2000)-1000;
         random3=double(rand()%2000)-1000;
         temp=sqrt(pow(random,2)+pow(random2,2)+pow(random3,2));
         random/=temp;
         random2/=temp;
         random3/=temp;
         sphere_dist[0][i]=acos(random3);
         if(random>=0 && random2>=0)
            sphere_dist[1][i]=atan(random2/random);

         else if(random<0 && random2>=0)
            sphere_dist[1][i]=Pi+atan(random2/random);

         else if(random<0 && random2<0)
            sphere_dist[1][i]=Pi+atan(random2/random);

         else if(random>0 && random2 <0)
            sphere_dist[1][i]=2*Pi+atan(random2/random);
      }
   }
   else
   {
      if(n_phi != 0)
         n_theta=n_points_sphere/n_phi;
      else
      {
         std::cout<<"CANNOT GENERATE A REGULAR SPHERICAL DISTRIBUTION WITH ZERO AZIMUTHAL ANGLE"<<std::endl;
         exit(EXIT_FAILURE);
      }
      for(int i=0;i!=n_theta;i++)
      {
         theta=i*acos(-1)/n_theta;
         for(int j=0;j!=n_phi;j++)
         {
            phi=j*2*acos(-1)/n_phi;
            sphere_dist[0][i*n_phi+j]=theta;
            sphere_dist[1][i*n_phi+j]=phi;
         }
      }
   }
}
