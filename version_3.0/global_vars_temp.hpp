#ifndef GLOBAL_VAR_H
#define GLOBAL_VAR_H

//Inlcude necessary headers
#include <cstdlib>
#include <string>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "mkl.h"
#include "omp.h"

//define physical and mathematical constants

const double au_to_fs(0.02418884);
const double pi(acos(-1));

//define the simulation parameters

const double total_time(25/au_to_fs);
const double h(0.001/au_to_fs);

//define the input file from which the electronic structure and the pulse parametersis infered
const std::string efield_path("e-field_input.in"); 

//define classes for electronic structure and electric field

class elec_struc  
{
   private:
      int m_n_states_neut;
      int m_n_states_cat;
      double *m_pot;
      double *m_dmx; 
      double *m_dmy; 
      double *m_dmz;

   public:
      elec_struc(std::string molpro_out_path,int n_states_neut,int n_states_cat);
      double dipole(int compos,int state_i,int state_f);
      double potential(int state);

};

class elec_field
{
   private:
      int m_number_of_pulse;
      double *m_pulse_energy;
      double *m_pulse_origin;
      double *m_pulse_intensity;
      double *m_pulse_sd;
      double *m_pulse_CEP;
      int    *m_pulse_pol;

   public:
      elec_field(std::string e_field_in_path);
      int pulse_number();
      double pulse_energy(int pulse_index);
      double pulse_origin(int pulse_index);
      double pulse_intensity(int pulse_index); 
      double pulse_sd(int pulse_index); 
      double pulse_CEP(int pulse_index);
      int pulse_pol(int pulse_index); 
      double efield_val(double vector[],double time);
};



#endif
