//
//  dyson_cube_writer.hpp
//  Dyson_computation
//
//  Created by Stephan van den Wildenberg on 18/01/17.
//  Copyright © 2017 Stephan van den Wildenberg. All rights reserved.
//

#ifndef dyson_cube_writer_hpp
#define dyson_cube_writer_hpp

#include <stdio.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

bool cube_header(int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int *n_nucl_dim,int *grid_size,int *num_of_nucl,int* basis_size,int *contraction_number,double nucl_coord,double **nucl_spher_pos,double *MO_coeff_neutral,double *dyson_mo_coeff,double **contraction_coeff,double **contraction_zeta,int* nucleus_basis_func,std::string *basis_func_type,int ** angular_mom_numbers,std::string Dyson_cube_loc,int state_neut,int state_cat,int n_nucl,int nx,int ny,int nz,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax);
bool cube_reader(int mo_index1,int mo_index2,int nx,int ny,int nz,std::string MO_cube_loc,double *cube_array);

#endif /* dyson_cube_writer_hpp */
