#include "photoion_comp.hpp"

int main(int argc,char* argv[])
{
   int photoion_comp(int argc, char* argv[]);
   omp_set_num_threads(8); 

   photoion_comp(argc,argv);

   return 0;
}

int photoion_comp(int argc, char* argv[])
{
    using namespace std;
    const bool symmetry(1);
    double temp_norm(0);
    const int n_sym(4);
    int n_states_neut(0);
    int n_states_neutral_sym[n_sym]={8,3,3,1};//{8,3,3,1};
    int n_states_cat(0);
    int n_states_cat_sym[n_sym]={2,1,1,0};//{2,1,1,0}
    int n_occ(0);
    int *n_occs;
    int n_closed(0);
    int *n_closeds;
    if(symmetry)
    {
       n_states_neut=0;
       n_states_cat=0;
       n_occs=new int [n_sym];
       n_closeds=new int [n_sym];
       for(int i=0;i!=n_sym;i++)
       {
          n_states_neut+=n_states_neutral_sym[i];
          n_states_cat+=n_states_cat_sym[i];
       }
    }
    std::cout<<n_states_neut<<" states in the neutral"<<std::endl<<n_states_cat<<" states in the cation"<<std::endl;
    int num_of_nucl(2);
    int basis_size(0);
    int* basis_size_sym=new int[n_sym];
    int **contraction_number;
    int ci_size_neut(0);
    int ci_size_cat(0);
    int ci_size_neut_sym[n_sym];
    int ci_size_cat_sym[n_sym];
    int n_elec_neut(4);//!!!! ecrire une routine qui cherche le nombre d'electrons dans l'output molpro!!!
   double Pi=acos(-1);
    double ***contraction_coeff;
    double ***contraction_zeta;
    int **contraction_number;
    int **nucl_basis_func;
    std::string** basis_func_type;


    //TEMPORARY VARIABLES
    double *overlap;
    double **mo_dipole;
    double *dyson_mo_basis_coeff;
    double *ci_vec_neut[2];
    double *ci_vec_cat[2];
    double sum(0);
    int index(0);
    int index2(0);
    int temp_int(0);

    string MO_cube_loc("MO_CUBE_LOC");
    string dyson_cube_loc("DYSON_CUBE_LOC");

//*****************************GET DATA FOR COMPUTING DYSON ORBITALS FROM MOLPRO OUTPUT FILE*****************************
        //GET THE NUMBER OF ELECTRONIC STATES AND THE SIZE OF THE ACTIVE SPACE
   if(symmetry)
   {
    size_query(n_occs,n_closeds,&basis_size, molpro_output_path,n_sym);
    n_occ=0;
    n_closed=0;
    for(int i=0;i!=n_sym;i++)
    {
       n_occ+=n_occs[i];
       n_closed+=n_closeds[i];
       std::cout<<"symmetry "<<i+1<<std::endl<<"closed "<<n_closeds[i]<<std::endl<<"occ "<<n_occs[i]<<std::endl;//DEBOGAGE 
    }
    std::cout<<"PROBE PI SIZE QUERY DONE"<<std::endl;//DEBOGAGE
   }
   else
      size_query(&n_occ,&n_closed,&basis_size,molpro_output_path);

    std::cout<<" number of occupied MO's : "<<n_occ<<std::endl;

    //GET THE SIZE OF THE CI VECTOR IN THE NEUTRAL AND THE CATION
    if(symmetry)
    {
    num_of_ci_reader(n_states_neutral_sym, n_states_cat_sym, ci_size_neut_sym, ci_size_cat_sym, molpro_output_path,n_occs,n_sym);
    ci_size_neut=0;
    ci_size_cat=0;
    for(int i=0;i!=n_sym;i++)
    {
       ci_size_neut+=ci_size_neut_sym[i];
       //std::cout<<ci_size_neut_sym[i]<<std::endl;
       ci_size_cat+=ci_size_cat_sym[i];
       //std::cout<<ci_size_cat_sym[i]<<std::endl;
    }
       std::cout<<"PROBE PI CI SIZE QUERY DONE"<<std::endl;//DEBOGAGE
    }
    else
      num_of_ci_reader(&n_states_neut, &n_states_cat, &ci_size_neut, &ci_size_cat, molpro_output_path,n_occs);
    std::cout<<"ci size neut is "<<ci_size_neut<<std::endl<<"ci size cat is "<<ci_size_cat<<std::endl;

    //ALLOCATE ARRAYS THAT DEPEND ON THE SIZE OF THE CI VECTOR
    ci_vec_neut[0]=new double[n_elec_neut*ci_size_neut+n_states_neut*ci_size_neut];//vector partitionned in two sections. Section 1 is filled with the mo label of each electron. section 2 is filled with CI coeff.
    ci_vec_neut[1]=new double[n_elec_neut*ci_size_neut];//this vector represents the spin state of each electron
    ci_vec_cat[0]=new double [(n_elec_neut-1)*ci_size_cat+n_states_cat*ci_size_cat];
    ci_vec_cat[1]=new double [(n_elec_neut-1)*ci_size_cat];
    dyson_mo_basis_coeff=new double[n_occ*n_states_neut*n_states_cat];

    for(int i=0;i!=n_elec_neut*ci_size_neut+n_states_neut*ci_size_neut;i++)
       ci_vec_neut[0][i]=0;
    for(int i=0;i!=(n_elec_neut-1)*ci_size_cat+n_states_cat*ci_size_cat;i++)
       ci_vec_cat[0][i]=0;

    std::cout<<"CI vector allocation done"<<std::endl;
    //GET THE CI COEFFICIENTS AND THE CONFIGURATIONS
    if(symmetry)
      ci_vec_reader(n_states_neutral_sym, n_states_cat_sym, n_occs,n_closeds, n_elec_neut,ci_size_neut_sym,ci_size_cat_sym, ci_vec_neut, ci_vec_cat, molpro_output_path,n_sym);
    else
      ci_vec_reader(&n_states_neut, &n_states_cat, &n_occ,&n_closed, n_elec_neut,&ci_size_neut,&ci_size_cat, ci_vec_neut, ci_vec_cat, molpro_output_path);
    std::cout<<"CI_VEC_READER ROUTINE ENDED WITHOUT ISSUE"<<std::endl;


//*****************************COMPUTE DYSON ORBITALS*****************************
    //COMPUTE THE OVERLAP MATRIX BETWEEN THE MO OF THE NEUTRAL AND THE MO OF THE CATION FROM THE AO OVERLAP MATRIX
    if(symmetry)
    overlap_MO(overlap,n_occs,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    else
    overlap_MO(overlap,&n_occ,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    //COMPUTE THE MOLECULAR ORBITALS TRANSITION DIPOLE MOMENT MATRIX
    if(symmetry)
    dipole_MO(mo_dipole,n_occs,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    else
    dipole_MO(mo_dipole,&n_occ,&basis_size,basis_size_sym,molpro_output_path,n_sym);

    //GET THE INFORMATION ABOUT THE BASIS SET
    contraction_coeff=new double**[n_sym];
    contraction_zeta=new double**[n_sym];
    contraction_number=new int*[n_sym];
    nucl_basis_func=new int*[n_sym];
    basis_func_type=new std::string*[n_sym];

    for(int s=0;s!=n_sym;s++)
    {
       contraction_number[s]=new int[basis_size_sym[s]];
       contraction_coeff[s]=new double*[basis_size_sym[s]];
       contraction_zeta[s]=new double*[basis_size_sym[s]];
       nucl_basis_func[s]=new int[basis_size_sym[s]];
       basis_func_type[s]=new std::string[basis_size_sym[s]]

    }
    basis_size_data_reader(n_sym, basis_size_sym,contraction_number,molpro_output_path);
    for(int s=0;s!=n_sym;s++)
    {
       for(int t=0;t!=basis_size_sym[s];t++)
       {
          contraction_coeff[s][t]=new double [contraction_number[s][t]];
          contraction_zeta[s][t]=new double [contraction_number[s][t]];
       }
    }
    basis_data_reader(n_sym,basis_size_sym,contraction_number,contraction_coeff,contraction_zeta,nucl_basis_func,basis_func_type,molpro_output_path);

    //COMPUTE THE DYSON MO COEFFICIENTS IN THE BASIS OF THE MO OF THE NEUTRAL
    dyson_mo_coeff_comp( n_states_neut,n_states_cat, n_occ,ci_size_neut, ci_size_cat, n_elec_neut, ci_vec_neut, ci_vec_cat,overlap, dyson_mo_basis_coeff);
    std::cout<<"DYSON COEFF ROUTINE ENDED WITHOUT ISSUE"<<std::endl;

    double dtemp(0);
    for(int a=0;a!=n_states_neut;a++)
    {
       dtemp=0;
       for(int k=0;k!=n_occ;k++)
       {
          dtemp+=dyson_mo_basis_coeff[a*n_states_cat*n_occ+0*n_occ+k]*dyson_mo_basis_coeff[a*n_states_cat*n_occ+0*n_occ+k];
       }
       std::cout<<"Norm of Dyson orbital "<<a<<"_"<<0<<" = "<<dtemp<<std::endl;
//cube_header(dyson_mo_basis_coeff,n_occ,n_states_neut,n_states_cat,neut_mo_cube_array,dyson_cube_loc.c_str(),a,0,2,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,mo_cube_array);
    }

    //BUILD THE CUBE OF THE DYSON ORBITALS
    //
}
