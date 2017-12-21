#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <complex>

double total_cs(int n_theta,int n_phi,double *theta,double *phi,double *mfpad);
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double prefactor,double *mfpad,double *redipole_x,double *imdipole_x,double *redipole_y,double *imdipole_y,double *redipole_z,double *imdipole_z);
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address);
bool angularly_resolved_dipole_reader(double *prefactor,double *theta,double *phi, double * redipole_x,double * imdipole_x,double * redipole_y,double * imdipole_y,double * redipole_z,double * imdipole_z,int n_points,std::string dipole_address);

int main(int argc, char *argv [])
{
   using namespace std;
   string dipole_address("/data1/home/stephan/LiH_fronuc_dyn/LiH_PICE_1.6125_7_0.txt");
   string spectrum_address("/data1/home/stephan/LiH_fronuc_dyn/LiH_CS_0_0.txt");
   int n_theta(90);
   int n_phi(120);
   int n_points(128*150);//(n_theta*n_phi);
   double theta_elec(0);//(acos(-1)/2);
   double phi_elec(0);//(acos(-1)/2);
   double x_comp(sin(theta_elec)*cos(phi_elec));
   double y_comp(sin(theta_elec)*sin(phi_elec));
   double z_comp(cos(theta_elec));
   double *theta=new double[n_points];
   double *phi=new double[n_points];
   double *redipole_x=new double[n_points];
   double *imdipole_x=new double[n_points];
   double *redipole_y=new double[n_points];
   double *imdipole_y=new double[n_points];
   double *redipole_z=new double[n_points];
   double *imdipole_z=new double[n_points];
   double *mfpad=new double[n_points];
   const double Pi(acos(-1));

   std::cout<<"Electric field components : (X-comp) , (Y-comp) , (Z-comp)"<<std::endl<<x_comp<<" , "<<y_comp<<", "<<z_comp<<std::endl;


   //angularly_resolved_dipole_reader(&prefactor,theta,phi,redipole_x,imdipole_x,redipole_y,imdipole_y,redipole_z,imdipole_z,n_points,dipole_address);

   /*for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
         x_comp=sin(theta[i*n_phi+j])*cos(phi[i*n_phi+j]);//X-component of the electric field
         y_comp=sin(theta[i*n_phi+j])*sin(phi[i*n_phi+j]);//Y-component of the electrif field
         z_comp=cos(theta[i*n_phi+j]);//Z-component of the electric field
        compute_mfpad(x_comp,y_comp,z_comp,n_points,prefactor,mfpad,redipole_x,imdipole_x,redipole_y,imdipole_y,redipole_z,imdipole_z);//Compute the angularly resolved differential cross section for the given direction of the electric field
        output<<theta[i*n_phi+j]<<"   "<<phi[i*n_phi+j]<<"   "<<total_cs(n_theta,n_phi,theta,phi,mfpad)<<std::endl;//Compute the total cross section for the given direction of the electric field
      }output<<std::endl;
   }*/
        //compute_mfpad(x_comp,y_comp,z_comp,n_points,prefactor,mfpad,redipole_x,imdipole_x,redipole_y,imdipole_y,redipole_z,imdipole_z);//Compute the angularly resolved differential cross section for the given direction of the electric field
        //output<<theta[i*n_phi+j]<<"   "<<phi[i*n_phi+j]<<"   "<<total_cs(n_theta,n_phi,theta,phi,mfpad)<<std::endl;//Compute the total cross section for the given direction of the electric field
         compute_spectrum( x_comp, y_comp, z_comp, n_points, dipole_address, spectrum_address);

   return 0;
}

bool angularly_resolved_dipole_reader(double *prefactor,double *theta,double *phi, double * redipole_x,double * imdipole_x,double * redipole_y,double * imdipole_y,double * redipole_z,double * imdipole_z,int n_points,std::string dipole_address)
{
   using namespace std;
   string tmp_str;

   ifstream dipole_file;
   dipole_file.open(dipole_address.c_str());
   if(!dipole_file.is_open())
   {
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   for(int i=0;i!=n_points;i++)
   {
      dipole_file>>tmp_str;
      dipole_file>>theta[i];
      dipole_file>>phi[i];
      if(i==0)
         dipole_file>>*prefactor;
      else
         dipole_file>>tmp_str;
      dipole_file>>redipole_x[i];
      dipole_file>>imdipole_x[i];
      dipole_file>>redipole_y[i];
      dipole_file>>imdipole_y[i];
      dipole_file>>redipole_z[i];
      dipole_file>>imdipole_z[i];
   }
   dipole_file.close();
   return 1;
   
}
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double prefactor,double *mfpad,double *redipole_x,double *imdipole_x,double *redipole_y,double *imdipole_y,double *redipole_z,double *imdipole_z)
{
   for(int i=0;i!=n_points;i++)
   {
      mfpad[i]=prefactor*(pow(x_comp*redipole_x[i]+y_comp*redipole_y[i]+z_comp*redipole_z[i],2)+pow(x_comp*imdipole_x[i]+y_comp*imdipole_y[i]+z_comp*imdipole_z[i],2));
   }
   return 1;
}
double total_cs(int n_theta,int n_phi,double *theta,double *phi,double *mfpad)
{
   double sum(0);
   double Pi=acos(-1);
   for(int i=0;i!=n_theta;i++)
   {
      for(int j=0;j!=n_phi;j++)
      {
          sum+=sin(theta[i*n_phi+j])*(Pi/n_theta)*(2*Pi/n_phi)*mfpad[i*n_phi+j];
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
   const int nk(150);
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
            cs[i]+=(4*pi/n_points)*(4*pi*pi*k[i*n_points/nk+j]/(137.0396096))*real(std::complex<double>(x_comp*RePICE_x[i*n_points/nk+j]+y_comp*RePICE_y[i*n_points/nk+j]+z_comp*RePICE_z[i*n_points/nk+j],x_comp*ImPICE_x[i*n_points/nk+j]+y_comp*ImPICE_y[i*n_points/nk+j]+z_comp*ImPICE_z[i*n_points/nk+j])*std::complex<double>(x_comp*RePICE_x[i*n_points/nk+j]+y_comp*RePICE_y[i*n_points/nk+j]+z_comp*RePICE_z[i*n_points/nk+j],-(x_comp*ImPICE_x[i*n_points/nk+j]+y_comp*ImPICE_y[i*n_points/nk+j]+z_comp*ImPICE_z[i*n_points/nk+j])));
         }
         std::cout<<pow(k[i*n_points/nk],2)*27.211/2<<", "<<cs[i]<<std::endl;
      }
      input.close(); 
   }
   return 1;
}
