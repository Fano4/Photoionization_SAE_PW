//
//  main.h
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 19/12/16.
//  Copyright © 2016 Stephan van den Wildenberg. All rights reserved.
//

#ifndef main_h
#define main_h

#include <iostream>
#include <string>
#include <iomanip>


#include "molpro_out_reader.hpp"
#include "dyson_cube_writer.hpp"
#include "computation.hpp"

int size_query(int* n_occ,int* basis_size,std::string molpro_out_path);
int overlap_MO(double matrix[],int* n_occ,int* basis_size,std::string molpro_out_path);
bool search(int *match_loc,std::string file_address,int research_from_this_pos,std::string pattern1, int num_of_entry_between_patterns12=0,std::string pattern2="",int num_of_entry_between_patterns23=0,std::string pattern3="");
int n_states_reader(int *n_states_neut,int *n_states_cat,std::string file_address);
int ci_vec_reader(int n_states_neut,int n_states_cat,int n_occ,int n_elec_neut,int ci_size_neut,int ci_size_cat,double **ci_vector_neut,double **ci_vector_cat,std::string file_address);
int num_of_ci_reader(int n_states_neut,int n_states_cat,int *n_ci_neut,int *n_ci_cat,std::string file_address);



#include "molpro_out_reader.cpp"
#include "dyson_cube_writer.cpp"
#include "computation.cpp"
#endif /* main_h */
