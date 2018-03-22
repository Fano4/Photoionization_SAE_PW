#include "global_vars.hpp"

elec_struc::elec_struc(std::string molpro_out_path,int n_states_neut,int n_states_cat)
{
   using namespace std;
    
   this->m_n_states_neut=n_states_neut;
   this->m_n_states_cat=n_states_cat;
    this->m_pot=new double[n_states_neut+n_states_cat];
    this->m_dmx=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];
    this->m_dmy=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];
    this->m_dmz=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];

    double *pot=new double[n_states_neut+n_states_cat];
    double *dmx=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];
    double *dmy=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];
    double *dmz=new double[(n_states_neut+n_states_cat)*(n_states_neut+n_states_cat)];   
    
    int curr_pos(0);
    curr_pos=potential_and_dipole_reader(molpro_out_path,0,n_states_neut,pot,dmx,dmy,dmz);
        if(curr_pos==0)
        {
           std::cout<<"PROBLEM WHILE READING POTENTIAL AND DIPOLE IN "<<molpro_out_path<<std::endl;
        }
        else
        {
            for (int k=0; k!=n_states_neut; k++)
            {   
                this->m_pot[k]=pot[k];
            }
            for(int k=0;k!=n_states_neut;k++)
            {
               for (int l=0; l!=n_states_neut;l++)
                {   
                   this->m_dmx[(n_states_neut+n_states_cat)*k+l]=dmx[(n_states_neut+n_states_cat)*k+l];
                   this->m_dmy[(n_states_neut+n_states_cat)*k+l]=dmy[(n_states_neut+n_states_cat)*k+l];
                   this->m_dmz[(n_states_neut+n_states_cat)*k+l]=dmz[(n_states_neut+n_states_cat)*k+l];
                }
            }
        }
    curr_pos=potential_and_dipole_reader(molpro_out_path,curr_pos,n_states_neut,pot,dmx,dmy,dmz);
        if(curr_pos==0)
        {
           std::cout<<"PROBLEM WHILE READING POTENTIAL AND DIPOLE IN "<<molpro_out_path<<std::endl;
        }
        else
        {
            for (int k=n_states_neut; k!=n_states_cat; k++)
            {   
                this->m_pot[k]=pot[k];
            }
            for(int k=n_states_neut;k!=n_states_neut+n_states_cat;k++)
            {
               for (int l=n_states_neut; l!=n_states_neut+n_states_cat;l++)
                {   
                   this->m_dmx[(n_states_neut+n_states_cat)*k+l]=dmx[(n_states_neut+n_states_cat)*k+l];
                   this->m_dmy[(n_states_neut+n_states_cat)*k+l]=dmy[(n_states_neut+n_states_cat)*k+l];
                   this->m_dmz[(n_states_neut+n_states_cat)*k+l]=dmz[(n_states_neut+n_states_cat)*k+l];
                }  
            } 
        }
}
double elec_struc::dipole(int compos,int state_i,int state_f)
{
    int temp(0);
    
    if(state_i>state_f)
    {
        temp=state_i;
        state_i=state_f;
        state_f=temp;
    }
   
    if((state_i<m_n_states_neut && state_f < m_n_states_neut)||(state_i > m_n_states_neut && state_f > m_n_states_neut))
    {
       switch (compos)
       {
           case 0:
               //name_indenter<<"dmx";
               return this->m_dmx[(m_n_states_neut+m_n_states_cat)*state_i+state_f];
               break;
           case 1:
               //std::cout<<this->m_dmy[n_states*state_i+state_f][x_index*gsize_y+y_index];
               return this->m_dmy[(m_n_states_neut+m_n_states_cat)*state_i+state_f];
               break;
           case 2:
               return this->m_dmz[(m_n_states_neut+m_n_states_cat)*state_i+state_f];
               //name_indenter<<"dmz";
               break;
       }
    }
    else
    {
    }
    return 0;
}


double elec_struc::potential(int state)
{
    return this->m_pot[state];
}

//Functions for reading molpro output file

elec_field::elec_field(std::string e_field_in_path)
{
   int n_pulses(0);
   double *energy;
   double *intensity;
   double *origin;
   double *sigma;
   double *CEP;
   int *ppol;

   if(!efield_param_reader( e_field_in_path, &n_pulses, &energy, &intensity, &origin, &sigma,&CEP,&ppol))
   {
      std::cout<<"FATAL ISSUE WHILE READING ELECTRIC FIELD PARAMETERS INPUT FILE "<<std::endl<<e_field_in_path<<std::endl;
      std::cout<<"CONTINUING WITHOUT ELECTRIC FIELD !!!!!"<<std::endl;
   }
   else
   {
      this->m_number_of_pulse=n_pulses;
      this->m_pulse_energy = new double[n_pulses];
      this->m_pulse_origin = new double[n_pulses];
      this->m_pulse_intensity = new double[n_pulses];
      this->m_pulse_sd = new double[n_pulses];
      this->m_pulse_CEP = new double[n_pulses];
      this->m_pulse_pol = new int[n_pulses];

      for(int i=0;i!=n_pulses;i++)
      {
         this->m_pulse_energy[i]=energy[i];
         this->m_pulse_intensity[i]=intensity[i];
         this->m_pulse_origin[i]=origin[i];
         this->m_pulse_sd[i]=sigma[i];
         this->m_pulse_CEP[i]=CEP[i];
         this->m_pulse_pol[i]=ppol[i];
      }
   }

}

double elec_field::pulse_energy(int pulse_index)
{
  return this->m_pulse_energy[pulse_index]; 
}
double elec_field::pulse_origin(int pulse_index)
{
   return this->m_pulse_origin[pulse_index];
}
double elec_field::pulse_intensity(int pulse_index)
{
   return this->m_pulse_intensity[pulse_index];
}
double elec_field::pulse_sd(int pulse_index)
{
   return this->m_pulse_sd[pulse_index];
}
double elec_field::pulse_CEP(int pulse_index)
{
   return this->m_pulse_CEP[pulse_index];
}
int elec_field::pulse_pol(int pulse_index)
{
   return this->m_pulse_pol[pulse_index];
}
double elec_field::efield_val(double vector[],double time)
{
    double E(0);
    double I(0);
    double t0(0);
    double S(0);
    double CEP(0);
     
    vector[0]=0;
    vector[1]=0;
    vector[2]=0;
    for(int i=0;i!=this->m_number_of_pulse;i++)
    {
      E=this->m_pulse_energy[i];
      I=this->m_pulse_intensity[i];
      t0=this->m_pulse_origin[i];
      S=this->m_pulse_sd[i];
      CEP=this->m_pulse_CEP[i];

      switch(m_pulse_pol[i])
      {
        case 0:
        vector[0]=I*exp(-pow((time-t0),2)/(2*S*S))*(cos(E*(time-t0)+CEP)-(time-t0)*sin(E*(time-t0)+CEP)/(E*S*S));;//X Component
             break;
        case 1:                                                                                                        
        vector[1]=I*exp(-pow((time-t0),2)/(2*S*S))*(cos(E*(time-t0)+CEP)-(time-t0)*sin(E*(time-t0)+CEP)/(E*S*S));//Y Component
             break;
        case 2:                                                                                                             
        vector[2]=I*exp(-pow((time-t0),2)/(2*S*S))*(cos(E*(time-t0)+CEP)-(time-t0)*sin(E*(time-t0)+CEP)/(E*S*S));//Z Component
             break;
      }
    }
    return 0;
}









