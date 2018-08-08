#ifndef hf5_photoion_hpp
#define hf5_photoion_hpp

#include <hdf5.h>
#include <H5Cpp.h>

bool write_output(std::string h5filename,int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int *n_nucl_dim,int *grid_size,int *num_of_nucl,int* basis_size,int *contraction_number,double *nucl_coord,double ***nucl_spher_pos,double ***mo_dipoles_mat,double **MO_coeff_neutral,double **dyson_mo_coeff,double **contraction_coeff,double **contraction_zeta,int* nucleus_basis_func,std::string *basis_func_type);

bool read_output(std::string h5filename,int* n_states_neut,int* n_states_cat,int *n_occ,int *n_closed,int *n_nucl_dim,int *grid_size,int *num_of_nucl,int* basis_size,int *contraction_number=NULL,double *nucl_coord=NULL,double ***nucl_spher_pos=NULL,double ***mo_dipoles_mat=NULL,double **MO_coeff_neutral=NULL,double **dyson_mo_coeff=NULL,double **contraction_coeff=NULL,double **contraction_zeta=NULL,int* nucleus_basis_func=NULL,std::string *basis_func_type=NULL);
int spherical_harmonics_translator(std::string basis_func_type,bool component);
std::string inverse_spherical_harmonics_translator(int l,int m);
   /*
    * The HDF5 PICE computation file will contain several groups and datasets that are listed below, respecting the groups hierarchy.
    *
    * ROOT
    * |
    * |--electronic_struc_parameters
    * |  |
    * |  |--symmetry_num_dataset
    * |  |
    * |  |--n_states_neut_dataset
    * |  |
    * |  |--n_states_cat_dataset
    * |  |
    * |  |--n_closed_tot_dataset
    * |  | 
    * |  |--n_occ_tot_dataset
    * |
    * |--nuclear_coordinates
    * |  |
    * |  |--n_nuclear_dim_dataset
    * |  |
    * |  |--coordinates_array_dataset
    * |  |  (n_nucl_dim X n_geom)
    * |  |
    * |  |--num_of_nucl_dataset
    * |  |
    * |  |--coordinates_atoms_cartesian_dataset
    * |  |  (n_geom X n_of_nucl X 3)
    * |
    * |--lcao_coeff_group
    * |  |
    * |  |--dyson_mo_coeff_matrix_dataset
    * |  |  (n_states_neut X n_states_cat X n_occ_tot X n_geom)
    * |  |
    * |  |--mo_coeff_neutral_dataset
    * |  |  ( n_occ_tot X basis_size X n_geom )
    * |
    * |--mo_dipole_group
    * |  |
    * |  |--mo_dipole_x
    * |  |  ((n_occ_tot+n_closed_tot) X (n_occ_tot+n_closed_tot) X n_geom)
    * |  |
    * |  |--mo_dipole_y
    * |  |  ((n_occ_tot+n_closed_tot) X (n_occ_tot+n_closed_tot) X n_geom)
    * |  |
    * |  |--mo_dipole_z
    * |  |  ((n_occ_tot+n_closed_tot) X (n_occ_tot+n_closed_tot) X n_geom)
    * |
    * |--basis_set_info_group 
    * |  |
    * |  |--basis_size_dataset
    * |  |
    * |  |--contraction_number_sym_dataset
    * |  |  (n_sym X basis_size)
    * |  |
    * |  |--contraction_coeff_sym_group
    * |  |  (there is n_sym datasets)
    * |  |  |
    * |  |  |--contraction_coeff_sym
    * |  |  |  (basis_size_sym X contraction_number_sym)
    * |  |
    * |  |--contraction_zeta_sym_group
    * |  |  (there is n_sym datasets)
    * |  |  |
    * |  |  |--contraction_zeta_sym
    * |  |  |  (basis_size_sym X contraction_number_sym)
    * |  |
    * |  |--nucl_basis_func_sym_group
    * |  |  (there is n_sym datasets)
    * |  |  |
    * |  |  |--nucl_basis_func_sym
    * |  |  |  (basis_size_sym)
    * |  |
    * |  |--basis_func_type_sym_group
    * |  |  (there is n_sym datasets)
    * |  |  |
    * |  |  |--basis_func_type_sym
    * |  |  |  (basis_size_sym X 2)
    * |  |  |  (for each basis functionm the spherical harmonics is given by its l and m_l numbers)
    * |  | 
    * |  |
    */

#endif /*hf5_photoion_hpp*/
