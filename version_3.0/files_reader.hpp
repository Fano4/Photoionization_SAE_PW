#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>


int n_states_reader(int *n_states_neut,int *n_states_cat,int *n_elec_neut,std::string file_address);
bool search(int *match_loc,std::string file_address,int research_from,std::string pattern1, int num_of_entry_between_patterns12,std::string pattern2,int num_of_entry_between_patterns23,std::string pattern3);
int ci_vec_reader(int *n_states_neut_s,int *n_states_cat_s,int *n_occ,int *n_closed,int n_elec_neut,int *ci_size_neut,int *ci_size_cat,double **ci_vector_neut,double **ci_vector_cat,std::string file_address,int n_sym);
int num_of_ci_reader(int *n_states_neut,int *n_states_cat,int *n_ci_neut,int *n_ci_cat,std::string file_address,int *n_occ,int n_sym);
int overlap_MO(double matrix[],int* n_occ,int* basis_size,int* basis_size_int,std::string molpro_out_path,int n_sym);
int size_query(int* n_occ,int *n_closed,int* basis_size,std::string molpro_out_path,int n_sym);
