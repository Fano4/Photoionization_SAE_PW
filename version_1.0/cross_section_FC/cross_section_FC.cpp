#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstring>

int fEPICE( double * PICE_energy,int n,double E);
void add_cs(int ini_index,int n,double k,double * theta,double *phi,std::complex<double> ** PICE,double factor,int comp_elec,double* value);

int main(int argc, char * argv[])
{
   using namespace std;
   string PICE_1_loc("HCN_PICE_0_0_test.txt");
   string PICE_2_loc("HCN_PICE_0_1_test.txt");
   string FC_1_loc("cationAp-Abs.vX001");
   string FC_2_loc("cationAs-Abs.vX001");
   string spectrum_loc("total_PICE_test.txt");
   string cs_spec_loc("sum_cs.txt");

   int length_FC_1(103);//number of lines in the FC input file
   int length_FC_2(37);
   int nk_PICE(100);
   int n_points_sphere(128);
   int length_PICE(nk_PICE*n_points_sphere);

   int spectrum_length(1000);
   double Emin(0.440998);
   double Emax(1.837492191);

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

   complex<double> ** total_PICE=new complex<double> * [3];
   total_PICE[0]=new complex<double> [spectrum_length*n_points_sphere];
   total_PICE[1]=new complex<double> [spectrum_length*n_points_sphere];
   total_PICE[2]=new complex<double> [spectrum_length*n_points_sphere];

   double *total_cs_x=new double[spectrum_length];
   double *total_cs_y=new double[spectrum_length];
   double *total_cs_z=new double[spectrum_length];

   double temp;
   double temp2;
   double E(0);

   ifstream input;
   ofstream output;
   ofstream output_cs;
   output.open(spectrum_loc.c_str());
   output.close();
   output_cs.open(cs_spec_loc.c_str());
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

   for(int e=0;e!=spectrum_length;e++)
   {
      E=Emin+e*(Emax-Emin)/spectrum_length;
      cout<<"index "<<e<<"/"<<spectrum_length<<" Energy "<<E<<endl;

      for(int i=0;i!=length_FC_1;i++)
      {
         if(E>=FC_1_energy[i])
         {
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i]), n_points_sphere,sqrt(2*(E-FC_1_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],PICE_1,FC_1_value[i],0,&total_cs_x[e]);
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i]), n_points_sphere,sqrt(2*(E-FC_1_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],PICE_1,FC_1_value[i],1,&total_cs_y[e]);
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i]), n_points_sphere,sqrt(2*(E-FC_1_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])],PICE_1,FC_1_value[i],2,&total_cs_z[e]);
            for(int j=0;j!=n_points_sphere;j++)
            {
                total_PICE[0][e*n_points_sphere+j]+=PICE_1[0][fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])+j]*FC_1_value[i];
                total_PICE[1][e*n_points_sphere+j]+=PICE_1[1][fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])+j]*FC_1_value[i];
                total_PICE[2][e*n_points_sphere+j]+=PICE_1[2][fEPICE(PICE_energy,length_PICE,E-FC_1_energy[i])+j]*FC_1_value[i];
            }
         }
      }
      for(int i=0;i!=length_FC_2;i++)
      {
         if(E>=FC_2_energy[i])
         {
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i]), n_points_sphere,sqrt(2*(E-FC_2_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],PICE_2,FC_2_value[i],0,&total_cs_x[e]);
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i]), n_points_sphere,sqrt(2*(E-FC_2_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],PICE_2,FC_2_value[i],1,&total_cs_y[e]);
                add_cs(fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i]), n_points_sphere,sqrt(2*(E-FC_2_energy[i])), &PICE_theta[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],&PICE_phi[fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])],PICE_2,FC_2_value[i],2,&total_cs_z[e]);
            for(int j=0;j!=n_points_sphere;j++)
            {
                total_PICE[0][e*n_points_sphere+j]+=PICE_2[0][fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])+j]*FC_2_value[i];
                total_PICE[1][e*n_points_sphere+j]+=PICE_2[1][fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])+j]*FC_2_value[i];
                total_PICE[2][e*n_points_sphere+j]+=PICE_2[2][fEPICE(PICE_energy,length_PICE,E-FC_2_energy[i])+j]*FC_2_value[i];
            }
         }
      }
      output.open(spectrum_loc.c_str(),ios_base::app);
      for (int j=0;j!=n_points_sphere;j++)
      {
          output<<sqrt(2*E)<<"  "<<PICE_theta[j]<<"  "<<PICE_phi[j]<<"  "<<real(total_PICE[0][e*n_points_sphere+j])<<"  "<<imag(total_PICE[0][e*n_points_sphere+j])<<"  "<<real(total_PICE[1][e*n_points_sphere+j])<<"  "<<imag(total_PICE[1][e*n_points_sphere+j])<<"  "<<real(total_PICE[2][e*n_points_sphere+j])<<"  "<<imag(total_PICE[2][e*n_points_sphere+j])<<endl;
      }output<<endl;
   output.close();
   output_cs.open(cs_spec_loc.c_str(),ios_base::app);
   output_cs<<E*27.211<<"  "<<total_cs_x[e]<<"  "<<total_cs_y[e]<<"  "<<total_cs_z[e]<<endl;
   output_cs.close();
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
