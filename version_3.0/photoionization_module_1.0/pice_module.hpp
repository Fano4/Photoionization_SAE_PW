#include "/data1/home/stephan/Photoionization_SAE_PW/version_3.0/algebra.cpp"
#include "/data1/home/stephan/Photoionization_SAE_PW/version_3.0/hf5_photoion.cpp"
#include "/data1/home/stephan/Photoionization_SAE_PW/version_3.0/Computation.cpp"

class pice_set {
   private:
       
      int *m_n_states_neut;
      int *m_n_states_cat;
      int *m_n_occ;
      int *m_n_closed;
      int *m_n_nucl_dim;
      int *m_grid_size;
      int *m_basis_size;
      int *m_num_of_nucl;
      int *m_max_cont_num;
      int* m_contraction_number;
      double **m_contraction_coeff;
      double **m_contraction_zeta;
      double **m_MO_coeff_neutral;
      double *m_nucl_coord;
      double ***m_nucl_spher_pos;
      double ***m_mo_dipole;
      double **m_dyson_mo_coeff;
      int* m_nucl_basis_func;
      std::string* m_basis_func_type;
/*
      std::complex<double> *k_part_s;
      std::complex<double> *k_part_p;
      std::complex<double> *k_part_d;
      std::complex<double> *k_part_f;
      std::complex<double> *k_part_s_gk;
      std::complex<double> *k_part_p_gk;
      std::complex<double> *k_part_d_gk;
      std::complex<double> *k_part_f_gk;
      std::complex<double> *k_part_s_gt;
      std::complex<double> *k_part_p_gt;
      std::complex<double> *k_part_d_gt;
      std::complex<double> *k_part_f_gt;
      std::complex<double> *k_part_s_gf;
      std::complex<double> *k_part_p_gf;
      std::complex<double> *k_part_d_gf;
      std::complex<double> *k_part_f_gf;
*/
   public:

      pice_set(std::string file_address);
      bool fill_pice(std::complex<double> *pice_x,std::complex<double> *pice_y,std::complex<double> *pice_z,int grid_index,int neut_st_index,int cat_st_index,double thet,double phi,double kp,double*pot_vec);
};
