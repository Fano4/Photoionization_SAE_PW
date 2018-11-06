#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <complex>

double total_cs(int n_theta,int n_phi,double *theta,double *phi,double *mfpad);
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double *mfpad,int kp,std::string dipole_address);
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address);
bool angularly_resolved_dipole_reader(double *prefactor,double *theta,double *phi, double * redipole_x,double * imdipole_x,double * redipole_y,double * imdipole_y,double * redipole_z,double * imdipole_z,int n_points,std::string dipole_address);

int main(int argc, char *argv [])
{
   using namespace std;
   string dipole_address("/data1/home/stephan/LiH_anion_cation_custom_6-311++G3df3dp/PICE_neut_0_0.txt");
   string spectrum_address("/data1/home/stephan/LiH_anion_cation_custom_6-311++G3df3dp/cs_anion_0_0.txt");
   int n_theta(20);
   int n_phi(40);
   int nk(50);
   int n_points(nk*2000);//nk*n_theta*n_phi);
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
//   int kp(5);
//        compute_mfpad(x_comp,y_comp,z_comp,n_points,mfpad,kp,dipole_address);//Compute the angularly resolved differential cross section for the given direction of the electric field
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
bool compute_mfpad(double x_comp,double y_comp,double z_comp,int n_points,double *mfpad,int kp,std::string dipole_address)
{
   int nk(50);
   int n_theta(20);
   int n_phi(40);
   double *RePICE_x=new double [n_points];
   double *ImPICE_x=new double [n_points];
   double *RePICE_y=new double [n_points];
   double *ImPICE_y=new double [n_points];
   double *RePICE_z=new double [n_points];
   double *ImPICE_z=new double [n_points];
   double *k=new double [n_points];
   double *thet=new double [n_points];
   double *phi=new double [n_points];

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
      for(int i=0;i!=nk*n_theta*n_phi;i++)
      {
         input>>k[i];
//         std::cout<<i<<"-"<<i%(n_theta*n_phi)<<"---"<<k[i]<<std::endl;
/*         if(i%(n_theta*n_phi)==0)
         {
            std::cout<<i/(n_theta*n_phi)<<"=>"<<k[i]<<std::endl;
         }*/
         input>>thet[i];
         input>>phi[i];
         input>>RePICE_x[i];
         input>>ImPICE_x[i];
         input>>RePICE_y[i];
         input>>ImPICE_y[i];
         input>>RePICE_z[i];
         input>>ImPICE_z[i];
      }
   input.close();
   for(int i=0;i!=n_theta*n_phi;i++)
   {
      mfpad[i]=real(std::complex<double>(x_comp*RePICE_x[kp*n_theta*n_phi+i]+y_comp*RePICE_y[kp*n_theta*n_phi+i]+z_comp*RePICE_z[kp*n_theta*n_phi+i],x_comp*ImPICE_x[kp*n_theta*n_phi+i]+y_comp*ImPICE_y[kp*n_theta*n_phi+i]+z_comp*ImPICE_z[kp*n_theta*n_phi+i])*std::complex<double>(x_comp*RePICE_x[kp*n_theta*n_phi+i]+y_comp*RePICE_y[kp*n_theta*n_phi+i]+z_comp*RePICE_z[kp*n_theta*n_phi+i],-(x_comp*ImPICE_x[kp*n_theta*n_phi+i]+y_comp*ImPICE_y[kp*n_theta*n_phi+i]+z_comp*ImPICE_z[kp*n_theta*n_phi+i])));
      std::cout<<k[kp*n_theta*n_phi+i]<<"   "<<thet[kp*n_theta*n_phi+i]<<"   "<<phi[kp*n_theta*n_phi+i]<<"   "<<mfpad[i]<<std::endl;
      if(i%n_phi==0 && i!=0)
         std::cout<<std::endl;
   }
   return 1;
}
}
bool compute_spectrum(double x_comp,double y_comp,double z_comp,int n_points,std::string dipole_address,std::string spectrum_address)
{
   int n_theta(20);
   int n_phi(40);
   int temp_int;
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
   const int nk(50);
   double total(0);
   double *cs=new double[nk];

   using namespace std;
   ifstream input;
   input.open("/data1/home/stephan/kxkykz.txt");
    double *distrib_cart=new double[3];

    for(int t=0;t!=2000;t++)
    {
       input>>temp_int;
       input>>distrib_cart[0];
       input>>distrib_cart[1];
       input>>distrib_cart[2];

       for(int p=0;p!=nk;p++)
       {
          theta[p*2000+t]=acos(distrib_cart[2]);
          phi[p*2000+t]=atan2(distrib_cart[1],distrib_cart[0]);
          if(phi[p*2000+t]<0)
             phi[p*2000+t]+=2*acos(-1);
       }
    }
    input.close();
    for(int p=0;p!=nk;p++)
    {
       for(int t=0;t!=2000;t++)
       {
          k[p*2000+t]=p*1.5/nk;
       }
    }

   input.open(dipole_address.c_str());

   if(!input.is_open())
   {      
      std::cout<<"CANNOT OPEN ANGULARLY RESOLVED IONIZATION DIPOLE FILE"<<std::endl<<"PROGRAM TERMINATION"<<std::endl;
      exit(EXIT_FAILURE);
   }
   else
   {
      for(int i=0;i!=nk*2000;i++)
      {
//         input>>k[i];
//         std::cout<<i<<"-"<<i%(n_theta*n_phi)<<"---"<<k[i]<<std::endl;
         if(i%(2000)==0)
         {
            std::cout<<i/(2000)<<"=>"<<k[i]<<std::endl;
         }
//         input>>theta[i];
//         input>>phi[i];
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
         for(int j=0;j!=2000;j++)
         {
//             for(int l=0;l!=n_phi;l++)
             {
//                std::cout<<k[i*n_theta*n_phi+j*n_phi+l]<<","<<theta[i*n_theta*n_phi+j*n_phi+l]<<","<<phi[i*n_theta*n_phi+j*n_phi+l]<<","<<std::endl;
      //          cs[i]+=pow(k[i*n_theta*n_phi+j*n_phi+l],2)*sin(theta[i*n_theta*n_phi+j*n_phi+l])*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)*real(std::complex<double>(x_comp*RePICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*RePICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*RePICE_z[i*n_theta*n_phi+j*n_phi+l],x_comp*ImPICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*ImPICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*ImPICE_z[i*n_theta*n_phi+j*n_phi+l])*std::complex<double>(x_comp*RePICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*RePICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*RePICE_z[i*n_theta*n_phi+j*n_phi+l],-(x_comp*ImPICE_x[i*n_theta*n_phi+j*n_phi+l]+y_comp*ImPICE_y[i*n_theta*n_phi+j*n_phi+l]+z_comp*ImPICE_z[i*n_theta*n_phi+j*n_phi+l])));
                cs[i]+=pow(k[i*2000+j],2)*(4*acos(-1)/2000)*real(std::complex<double>(x_comp*RePICE_x[i*2000+j]+y_comp*RePICE_y[i*2000+j]+z_comp*RePICE_z[i*2000+j],x_comp*ImPICE_x[i*2000+j]+y_comp*ImPICE_y[i*2000+j]+z_comp*ImPICE_z[i*2000+j])*std::complex<double>(x_comp*RePICE_x[i*2000+j]+y_comp*RePICE_y[i*2000+j]+z_comp*RePICE_z[i*2000+j],-(x_comp*ImPICE_x[i*2000+j]+y_comp*ImPICE_y[i*2000+j]+z_comp*ImPICE_z[i*2000+j])));
             }
         }
//         total+=cs[i]*(k[nk*n_theta*n_phi-1]-k[0])/nk;
         total+=cs[i]*(k[nk*2000-1]-k[0])/nk;
//         std::cout<<k[i*2000]<<", "<<cs[i]<<" => total = "<<total<<std::endl;
         std::cout<<k[i*2000]*k[i*2000]*27.211/2<<"   "<<cs[i]*1.152321368*pow(0.592e-10,2)*10000<<std::endl;
 //        std::cout<<pow(k[i*n_theta*n_phi],2)*27.211/2<<", "<<cs[i]<<std::endl;
         //std::cout<<pow(k[i*n_theta*n_phi],2)*27.211/2<<", "<<cs[i]*Lx*Ly*Lz/pow((2*pi),3)<<std::endl;
      }
      input.close(); 
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
