#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstring>

int fEPICE( double * PICE_energy,int n,double E);
void add_cs(int ini_index,int n,double k,double * theta,double *phi,std::complex<double> ** PICE,double factor,int comp_elec,double* value);

int main(int argc, char* argv[])
{
   using namespace std;


   const double total_energy(0.7790966888);
   const double uncertainty(0.001/27.211);
   const int pes_length(8000);
   const double Kmin(0);
   const double Kmax(total_energy);//The kinetic energy of the electron is bound by the energy of the photon
   double *pes_x=new double [pes_length];//the abcissa of the pes is the kinetic energy of the electron
   double *pes_y=new double [pes_length];//the abcissa of the pes is the kinetic energy of the electron
   double *pes_z=new double [pes_length];//the abcissa of the pes is the kinetic energy of the electron
   double K(0);

   double temp(0);
   double temp2(0);

   string PICE_1_loc("HNC_PICE_0_0_large.txt");
   string PICE_2_loc("HNC_PICE_0_1_large.txt");
   string FC_1_loc("cationAp-Abs.vX011");
   string FC_2_loc("cationAs-Abs.vX011");
   string pes_loc("HNC_photoelec_spectrum.txt");
   int length_FC_1(103);//number of lines in the FC input file
   int length_FC_2(37);
   int nk_PICE(1000);
   int n_points_sphere(256);
   int length_PICE(nk_PICE*n_points_sphere);

   double * FC_1_energy=new double[length_FC_1];
   double * FC_2_energy=new double[length_FC_2];
   double * FC_1_value=new double[length_FC_1];
   double * FC_2_value=new double[length_FC_2];
   double * PICE_energy=new double[length_PICE];
   double * PICE_theta=new double[length_PICE];
   double * PICE_phi=new double[length_PICE];

   complex<double> ** PICE_1=new complex<double> * [3];
   PICE_1[0]=new complex<double> [length_PICE];
   PICE_1[1]=new complex<double> [length_PICE];
   PICE_1[2]=new complex<double> [length_PICE];
   complex<double> ** PICE_2=new complex<double> * [3];
   PICE_2[0]=new complex<double> [length_PICE];
   PICE_2[1]=new complex<double> [length_PICE];
   PICE_2[2]=new complex<double> [length_PICE];


   ifstream input;
   ofstream output;
   output.open(pes_loc.c_str());
   output.close();

   input.open(FC_1_loc.c_str());
   if(!input.is_open())
   {
      std::cout<<"ERROR: CANNOT OPEN INPUT"<<endl<<FC_1_loc.c_str()<<endl<<"PLEASE CHECK FILE"<<endl; 
   }
   else
   {
      for(int i=0;i!=length_FC_1;i++)
      {
         input>>FC_1_energy[i];
         input>>temp;
         input>>temp;
         input>>temp;
         input>>temp;
         input>>temp;
         FC_1_value[i]=temp*temp;
         input>>temp;
      }
         input.close();
   }
   //for(int i=0;i!=length_FC_1;i++)
  // {
  //    std::cout<<FC_1_energy[i]<<"  "<<FC_1_value[i]<<endl;
  // }
//return 0;
   input.open(FC_2_loc.c_str());
   if(!input.is_open())
   {
      std::cout<<"ERROR: CANNOT OPEN INPUT"<<endl<<FC_2_loc.c_str()<<endl<<"PLEASE CHECK FILE"<<endl; 
   }
   else
   {
      for(int i=0;i!=length_FC_2;i++)
      {
         input>>FC_2_energy[i];
         input>>temp;
         input>>temp;
         input>>temp;
         input>>temp;
         input>>temp;
         FC_2_value[i]=temp*temp;
         input>>temp;
      }
         input.close();
   }
//   for(int i=0;i!=length_FC_2;i++)
//   {
//      std::cout<<FC_2_energy[i]<<"  "<<FC_2_value[i]<<endl;
//   }
//return 0;
   input.open(PICE_1_loc.c_str());
   if(!input.is_open())
   {
      std::cout<<"ERROR: CANNOT OPEN INPUT"<<endl<<PICE_1_loc.c_str()<<endl<<"PLEASE CHECK FILE"<<endl; 
   }
   else
   {
      for(int i=0;i!=length_PICE;i++)
      {
         input>>temp;
         PICE_energy[i]=pow(temp,2)/2;
         input>>PICE_theta[i];
         input>>PICE_phi[i];
         input>>temp;
         input>>temp2;
         PICE_1[0][i]=complex<double>(temp,temp2);
         input>>temp;
         input>>temp2;
         PICE_1[1][i]=complex<double>(temp,temp2);
         input>>temp;
         input>>temp2;
         PICE_1[2][i]=complex<double>(temp,temp2);
      }
         input.close();
   }

   input.open(PICE_2_loc.c_str());
   if(!input.is_open())
   {
      std::cout<<"ERROR: CANNOT OPEN INPUT"<<endl<<PICE_2_loc.c_str()<<endl<<"PLEASE CHECK FILE"<<endl; 
   }
   else
   {
      for(int i=0;i!=length_PICE;i++)
      {
         input>>temp;
         input>>temp;
         input>>temp;
         
         input>>temp;
         input>>temp2;
         PICE_2[0][i]=complex<double>(temp,temp2);
         input>>temp;
         input>>temp2;
         PICE_2[1][i]=complex<double>(temp,temp2);
         input>>temp;
         input>>temp2;
         PICE_2[2][i]=complex<double>(temp,temp2);
      }
         input.close();
   }


   //we scan the kinetic energies of the electron from Kmin to Kmax.
   //K=0 kinetic energy means that the total energy was absorbed by the ionization and vibrational excitation
   //K=total_energy means that the ionization has absorbed no energy
   for(int k=0;k!=pes_length;k++)
   {
      K=Kmin+k*(Kmax-Kmin)/pes_length;
      //for each kinetic energy, we scan which ionization/vibrational excitation can emit an electron with this kinetic energy so that the energy is conserved.
      for(int e=0;e!=length_FC_1;e++)
      {
         if(FC_1_energy[e]+K<=total_energy+uncertainty &&FC_1_energy[e]+K>=total_energy-uncertainty)//when kinetic energy plus electronic/vibrational energy is equal to total energy, then this channel is open and the cross section is added to the spectrum
         {
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_1_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],PICE_1,FC_1_value[e],0,&pes_x[k]);
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_1_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],PICE_1,FC_1_value[e],1,&pes_y[k]);
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_1_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_1_energy[e])],PICE_1,FC_1_value[e],2,&pes_z[k]);
         }
         else
            continue;
      }
      for(int e=0;e!=length_FC_2;e++)
      {
         if(FC_2_energy[e]+K<=total_energy+uncertainty && FC_2_energy[e]+K>=total_energy-uncertainty)//when kinetic energy plus electronic/vibrational energy is equal to total energy, then this channel is open
         {
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_2_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],PICE_2,FC_2_value[e],0,&pes_x[k]);
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_2_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],PICE_2,FC_2_value[e],1,&pes_y[k]);
                add_cs(fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e]), n_points_sphere,sqrt(2*(total_energy-FC_2_energy[e])), &PICE_theta[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],&PICE_phi[fEPICE(PICE_energy,length_PICE,total_energy-FC_2_energy[e])],PICE_2,FC_2_value[e],2,&pes_z[k]);
         }
         else
            continue;
      }
      output.open(pes_loc.c_str(),ios_base::app);
      output<<(total_energy-K)*27.211<<"  "<<pes_x[k]<<"  "<<pes_y[k]<<"  "<<pes_z[k]<<endl;
      output.close();
   }

   return 0;
}
int fEPICE( double * PICE_energy,int n,double E)
{
   if(E<PICE_energy[0])
   {
      return 0;
   }
   if(E>PICE_energy[n-1])
   {
      return n-1;
   }
   for(int i=0;i!=n;i++)
   {
      if(E>=PICE_energy[i] && E<PICE_energy[i+1])
         return i;
      else 
         continue;
   }
   std::cout<<"ERROR : ENERGY NOT MAPPED ON TABLE"<<std::endl;
   return -1;
}
void add_cs(int ini_index,int n,double k,double * theta,double *phi,std::complex<double> ** PICE,double factor,int comp_elec,double *value)
{

   const double pi(acos(-1));
   const double Lx(26.6337);
   const double Ly(26.6337);
   const double Lz(27.876);
   double theta_elec;
   double phi_elec;
   if(comp_elec==0)
   {
      theta_elec=pi/2;
      phi_elec=0;
   }
   else if(comp_elec==1)
   {
      theta_elec=pi/2;
      phi_elec=pi/2;
   }
   else if(comp_elec==2)
   {
      theta_elec=0;
      phi_elec=0;
   }
   double x_comp(sin(theta_elec)*cos(phi_elec));
   double y_comp(sin(theta_elec)*sin(phi_elec));
   double z_comp(cos(theta_elec));

         for(int j=0;j!=n;j++)
         {
            *value+=factor*(4*pi/n)*(2*k*Lx*Ly*Lz/(137.0396096*pi))*real(std::complex<double>(x_comp*std::real(PICE[0][ini_index+j])+y_comp*std::real(PICE[1][ini_index+j])+z_comp*std::real(PICE[2][ini_index+j]),x_comp*std::imag(PICE[0][ini_index+j])+y_comp*std::imag(PICE[1][ini_index+j])+z_comp*std::imag(PICE[2][ini_index+j]))*std::complex<double>(x_comp*std::real(PICE[0][ini_index+j])+y_comp*std::real(PICE[1][ini_index+j])+z_comp*std::real(PICE[2][ini_index+j]),-(x_comp*std::imag(PICE[0][ini_index+j])+y_comp*std::imag(PICE[1][ini_index+j])+z_comp*std::imag(PICE[2][ini_index+j]))));
         }
}
