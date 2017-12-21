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
#include <fstream>
#include <sstream>
#include <iomanip>

#include "dyson_vars.hpp"
#include "molpro_out_reader.hpp"


bool cube_header(double *dyson_MO_basis_coeff,int n_occ,int n_states_neut,int n_states_cat,std::string MO_cube_loc,std::string Dyson_cube_loc);
bool cube_reader(int mo_index,int nx,int ny,int nz,std::string MO_cube_loc,double *cube_array);

#endif /* dyson_cube_writer_hpp */
