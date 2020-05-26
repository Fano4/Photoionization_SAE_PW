#include <vector>
double azim_integ(int m1);
double two_azim_integ(int m1,int m2);
double three_azim_integ(int m1,int m2,int m3);

void Jint_sort_indices(int* l1,int* l2,int* l3,int* m1,int* m2,int* m3);
double Jint_signflip_renormalize(int l1,int l2,int l3,int* m1,int* m2,int* m3);
double ALP_normalize(int l,int m);
double Jint_normalize(int l1,int l2,int l3,int m1,int m2,int m3);
bool Jint_special_cases(int l1,int l2,int l3,int m1,int m2,int m3,double* result);
double ALP_integral(int l,int m);
double two_ALP_integral(int l1,int l2,int m1,int m2);
double three_ALP_J_integral(int l1,int l2,int l3,int m1,int m2,int m3);
double three_Ylm_integ(int l1,int l2,int l3,int m1,int m2,int m3);


double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_p1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1_D(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_p1_D(int l1,int l2,int l3,int m1,int m2,int m3);

double I_m1_integral(int m1,int m2,int m3);
double I_p1_integral(int m1,int m2,int m3);
double I_m1_D_integral(int m1,int m2,int m3);
double I_p1_D_integral(int m1,int m2,int m3);

//double spherical_overlap_integral(double xa,double xb,double m,double p);
//double obara_saika_ovlp(double xa,double xb,double zeta_a,double zeta_b,int la,int lb);
//double bessel_gaussian_poly_integral(unsigned int l1,unsigned int l2,double m,double r);
double prim_radial_ovlp(unsigned int la,unsigned int lb,unsigned int l,double zet_a,double zet_b,double r);
double prim_ovlp(std::vector<double> ra,std::vector<double> rb,double zeta_a,double zeta_b,unsigned int la,unsigned int lb,int ma,int mb);
void ao_ovlp(std::vector<double> ra,std::vector<double> rb,std::vector<int> nuc_bas_fun_a,std::vector<int> nuc_bas_fun_b,std::vector<int> cont_num_a,std::vector<int> cont_num_b,std::vector<double> zet_a,std::vector<double> zet_b,std::vector<double> cont_coeff_a, std::vector<double> cont_coeff_b,std::vector<unsigned int> la,std::vector<unsigned int> lb,std::vector<int> ma,std::vector<int> mb,std::vector<double>* S);
void MO_ovlp(std::vector<double> S,std::vector<double> lcao_a,std::vector<double> lcao_b,std::vector<double>* MO_S);
double slater_ovlp(std::vector<int> mo_a,std::vector<int> mo_b,std::vector<int> spin_a,std::vector<int> spin_b,std::vector<double> MO_S);


double determinant(double *A,int dim);
void transpose(double *A,double *B, int dim1, int dim2);
void matrix_product(double *C,double *A,double *B,int dim1,int dim2,int dim3);
double prefactor_rYlm(int l,int m);
double rYlm (int l,int m,double thet,double phi);

void compute_bessel_pice_mo(double*** pice_ortho_mo,double*** pice_ddx_mo,double*** pice_ddy_mo,double*** pice_ddz_mo,int jl_max,int n_occ,int basis_size,int nk,double kmax,double *MO_coeff_neutral,double **contraction_zeta,double **contraction_coeff,int * contraction_number,double** nucl_spher_pos,int *nucl_basis_func,int** angular_mom_numbers);
