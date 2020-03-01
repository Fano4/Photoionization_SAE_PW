#include "photoion_comp.hpp"

#include <iomanip>

int main(int argc,char* argv[])
{
   int photoion_comp(int argc, char* argv[]);
   omp_set_num_threads(1); 
   int lmax=10;
   for(int l1=0;l1!=lmax+1;l1++)
   {
      for(int m1=0;m1!=l1+1;m1++)
      {
//         std::cout<<l1*l1+l1+m1<<"/120"<<std::endl;;
//         std::cout<<"{";
         for(int l2=0;l2!=lmax+1;l2++)
         {
            for(int m2=0;m2!=l2+1;m2++)
            {
//               std::cout<<"{";
               for(int l3=0;l3!=lmax+1;l3++)
               {
                  for(int m3=0;m3!=l3+1;m3++)
                  {
                     std::cout<<"<<<<"<<l1<<","<<l2<<","<<l3<<";"<<m1<<","<<m2<<","<<m3<<std::endl;

                     std::cout<<J_int_m2(l1,l2,l3,m1,m2,m3)<<" == "<<test2_integral(l1,l2,l3,m1,m2,m3)<<std::endl; // wigner3j(l1,l2,l3,m1,m2,m3);
                      /*  
                     if(l3*l3+l3+m3<lmax*lmax+2*lmax)
                     std::cout<<std::setprecision(15)<<wigner3j(l1,l2,l3,m1,m2,m3)<<",";
                     else
                     std::cout<<std::setprecision(15)<<wigner3j(l1,l2,l3,m1,m2,m3);
          * /           
//                     if(l3*l3+l3+m3<lmax*lmax+2*lmax)
                     if((l3*(l3+1)/2+m3)<(lmax*(lmax+1)/2+lmax))
                        std::cout<<std::setprecision(15)<<gaunt_formula(l1,l2,l3,m1,m2,m3)<<",";
                     else
                        std::cout<<std::setprecision(15)<<gaunt_formula(l1,l2,l3,m1,m2,m3);
                    */ 
                  }
               }
//               if(l2*l2+l2+m2<lmax*lmax+2*lmax)
//               if((l2*(l2+1)/2+m2)<lmax*(lmax+1)/2+lmax)
//                   std::cout<<"},"<<std::endl;
//               else
//                  std::cout<<"}";
            }
         }
//               if((l1*(l1+1)/2+m1)<lmax*(lmax+1)/2+lmax)
//                   std::cout<<"},"<<std::endl;
//               else
//                  std::cout<<"}";
      }
   }
   exit(EXIT_SUCCESS);
   double* r0=new double [3];

   r0[0]=0;
   r0[1]=0;
   r0[2]=0;

  /* 
   for(int l1=0;l1!=5;l1++)
   {
   for(int m1=0;m1!=l1+1;m1++)
   {
   for(int l2=0;l2!=5;l2++)
   {
   for(int m2=0;m2!=l2+1;m2++)
   {
   for(int l3=0;l3!=5;l3++)
   {
   for(int m3=0;m3!=l3+1;m3++)
   {
   std::cout<<l1<<","<<m1<<" - "<<l2<<","<<m2<<" - "<<l3<<","<<m3<<" : JP1D = "<<J_int_p1_D(l1,l2,l3,m1,m2,m3)<<std::endl;
   }
   }
   }
   }
   }
   }
   exit(EXIT_SUCCESS);
*/
//   pw_bessel_comparison(0.3,0.723,2.3,0.3430,1.145,4.2049);
//   pw_bessel_overlap_comparison(1,-1,0.125,0.3,0.5,2.1,r0);
   pw_bessel_gradient_y_comparison(1,-1,0.384,0.1,0.22,2.1,r0);
//   pw_bessel_comparison(0.3,0.25,2.1,3.75,1.25,0.887);
//   for(int k=0;k!=256;k++)
   {
//      std::cout<<(k+1)*15/2560.<<","<<j_l(1,(k+1)*15/2560.,lnfact_memo)<<std::endl;
//      numerical_integral(1,1,1,1,0.125,(k+1)*15/2560.,r0);
//      analytic_integral(1,1,1,1,0.125,(k+1)*15/2560.,r0);


   }
   exit(EXIT_SUCCESS);
   
   /*
   int l1(0);
   int l2(0);
   int l3(0);
   int m1(0);
   int m2(0);
   int m3(0);
   double val(0);
   double aval(0);
   double begin(0);
   double end(0);

   for(int l1=0;l1!=5;l1++)
   {
      for(int l2=0;l2!=4;l2++)
      {
         for(int l3=0;l3!=l1+l2+1;l3++)
         {
            for(int m1=0;m1!=l1+1;m1++)
            {
               for(int m2=0;m2!=l2+1;m2++)
               {
                  for(int m3=0;m3!=l3+1;m3++)
                  {
    begin=omp_get_wtime();
    val=(test2_integral(l1,l2,l3,m1,m2,m3,lnfact_memo));
    end=omp_get_wtime(); 
    std::cout<<std::endl<<std::fixed<<std::setprecision(6)<<(end-begin)/10<<"s is the time with a precision of"<<omp_get_wtick()<<std::endl;
    begin=omp_get_wtime();
    aval=(gaunt_formula(l1,l2,l3,m1,m2,m3,lnfact_memo));
    end=omp_get_wtime(); 
    std::cout<<std::endl<<std::fixed<<std::setprecision(6)<<(end-begin)/10<<"s is the time with a precision of"<<omp_get_wtick()<<std::endl;
   std::cout<<l1<<","<<l2<<","<<l3<<","<<m1<<","<<m2<<","<<m3<<"***"<<val<<std::endl;
   std::cout<<l1<<","<<l2<<","<<l3<<","<<m1<<","<<m2<<","<<m3<<"+++"<<aval<<std::endl;
                  }
               }
            }
         }
      }
   }

   exit(EXIT_SUCCESS);
*/
   photoion_comp(argc,argv);

   return 0;
}

int photoion_comp(int argc, char* argv[])
{
    using namespace std;
    //variables independent of grid size
    const bool symmetry(1);
    const int n_sym(4);
    int nucl_dim(1);
    int n_states_neut(0);
    int n_states_neutral_sym[n_sym]={10,4,4,1};//{8,3,3,1};
    int n_states_cat(0);
    int n_states_cat_sym[n_sym]={4,1,1,0};//{2,1,1,0}
    int n_occ(0);
    int *n_occs;
    int n_closed(0);
    int *n_closeds;
    int num_of_nucl(2);
    int basis_size(0);
    int* basis_size_sym=new int[n_sym];
    int ci_size_neut(0);
    int ci_size_cat(0);
    int ci_size_neut_sym[n_sym];
    int ci_size_cat_sym[n_sym];
    double **ci_vec_neut=new double*[2];
    double **ci_vec_cat=new double*[2];
    double ***contraction_coeff_sym;
    double ***contraction_zeta_sym;
    int **contraction_number_sym;
    int **nucl_basis_func_sym;
    std::string** basis_func_type_sym;
    int ***angular_mom_numbers_sym;
    int **angular_mom_numbers;
    double **contraction_coeff;
    double **contraction_zeta;
    double **contraction_coeff_array;
    double **contraction_zeta_array;
    int *contraction_number;
    int *nucl_basis_func;
    std::string* basis_func_type;
    int n_elec_neut(4);//!!!! ecrire une routine qui cherche le nombre d'electrons dans l'output molpro!!!
    int max_contraction_num(0);


    double testtime=omp_get_wtime();
    std::cout<<std::fixed<<std::setprecision(6)<<"test_time = "<<testtime<<"s is the time with a precision of"<<omp_get_wtick()<<std::endl;
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
    //variables depending on grid size
    const std::string hf5_outfile("LiH_testing_continuum.h5");
    int grid_size(1);
    double **MO_coeff_neutral=new double*[grid_size];
    double *overlap;
    double ***mo_dipole=new double**[grid_size];
    double ***mo_dipole_spher=new double**[grid_size];
    double **dyson_mo_basis_coeff=new double*[grid_size];
    double ***tran_den_mat_mo=new double ** [grid_size];
    double xp(0);
    double rLi(0);
    double rH(0);
    double *nucl_coord=new double[grid_size];

    double ***nucl_spher_pos=new double**[grid_size];

    //TEMPORARY VARIABLES
    double sum(0);
    int index(0);
    int index2(0);
    int temp_int(0);
    double xmin(1.6125);
    double xmax(1.6125);
    double mLi(6.015122795);
    double mH(1.007825);
    std::string file_root("/data1/home/stephan/LiH_gridtest_+++custom/LiH_");
    stringstream ss_molpro_file;
    std::string molpro_output_path;

//       std::cout<<" Factorial of "<<" = "<<exp(ln_factorial(25,lnfact_memo)-ln_factorial(20,lnfact_memo))<<std::endl;
    /*
     *
     *
     * 

    int pl1(2);
    int pl2(1);
    int pl3(1);
    int pm1(2);
    int pm2(1);
    int pm3(1);

    std::cout<<gaunt_formula(pl1+1,pl2,pl3,pm1+1,pm2,pm3,lnfact_memo)<<std::endl;
    std::cout<<gaunt_formula(pl1+1,pl2,pl3,pm1-1,pm2,pm3,lnfact_memo)<<std::endl;
    std::cout<<J_int_m2(pl1,pl2,pl3,pm1,pm2,pm3,lnfact_memo)<<std::endl;
    *
     *
     *
     * 
    exit(EXIT_SUCCESS);
*/
    for(int x=0;x!=grid_size;x++)
    {
       nucl_spher_pos[x]=new double* [num_of_nucl]; 


       for(int i=0;i!=num_of_nucl;i++)
       {
          nucl_spher_pos[x][i]=new double [3]; 
          for(int d=0;d!=3;d++)
          {
              nucl_spher_pos[x][i][d]=0;
          }
       }
    }

//*****************************GET DATA FOR COMPUTING DYSON ORBITALS FROM MOLPRO OUTPUT FILE*****************************
for(int x=0;x!=grid_size;x++)
{

   xp=xmin+x*(xmax-xmin)/grid_size;
   rH=mLi*xp/(mLi+mH);
   rLi=mH*xp/(mLi+mH);

   nucl_coord[x]=xp;
   ss_molpro_file.str("");
   ss_molpro_file<<file_root.c_str()<<xp<<".out";
   molpro_output_path=ss_molpro_file.str();

   std::cout<<"WORKING WITH MOLPRO FILE \""<<molpro_output_path.c_str()<<"\""<<std::endl;

   //*****************************THIS SECTION OF THE CODE DOES NOT DEPEND ON NUCLEAR POSITIONS. IT MUST BE CALLED ONLY ONCE*****************************
   //
   //
   if(x==0)
   {
      if(symmetry)
      {
         size_query(n_occs,n_closeds,basis_size_sym, molpro_output_path,n_sym);
         n_occ=0;
         n_closed=0;
         basis_size=0;
         for(int i=0;i!=n_sym;i++)
         {
            n_occ+=n_occs[i];
            n_closed+=n_closeds[i];
            basis_size+=basis_size_sym[i];
            std::cout<<"symmetry "<<i+1<<std::endl<<"closed "<<n_closeds[i]<<std::endl<<"occ "<<n_occs[i]<<std::endl;//DEBOGAGE 
         }
         std::cout<<"PROBE SIZE QUERY DONE"<<std::endl;//DEBOGAGE
      }
      else
         size_query(&n_occ,&n_closed,&basis_size,molpro_output_path);

   //   std::cout<<" number of occupied MO's : "<<n_occ<<std::endl;
    //GET THE INFORMATION ABOUT THE BASIS SET
    contraction_coeff_sym=new double**[n_sym];
    contraction_zeta_sym=new double**[n_sym];
    contraction_number_sym=new int*[n_sym];
    nucl_basis_func_sym=new int*[n_sym];
    basis_func_type_sym=new std::string*[n_sym];
    angular_mom_numbers_sym=new int**[n_sym];

    std::cout<<"INITIALIZED SYMMETRY DEPENDENT ARRAYS. INITIALIZING BASIS SET SIZE DEPENDENT ARRAYS :"<<std::endl;
    int total(0);

    for(int s=0;s!=n_sym;s++)
    {
       std::cout<<"size sym "<<s<<" = "<<basis_size_sym[s]<<std::endl;
       contraction_number_sym[s]=new int[basis_size_sym[s]];
       contraction_coeff_sym[s]=new double*[basis_size_sym[s]];
       contraction_zeta_sym[s]=new double*[basis_size_sym[s]];
       nucl_basis_func_sym[s]=new int[basis_size_sym[s]];
       basis_func_type_sym[s]=new std::string[basis_size_sym[s]];
       angular_mom_numbers_sym[s]=new int*[basis_size_sym[s]];
       total+=basis_size_sym[s];
    }
//    std::cout<<"probe"<<std::endl;
    contraction_coeff=new double*[total];
    contraction_zeta=new double*[total];
    contraction_number=new int[total];
    nucl_basis_func=new int[total];
    basis_func_type=new std::string[total];
    angular_mom_numbers=new int*[total];

    basis_size_data_reader(n_sym, basis_size_sym,contraction_number_sym,molpro_output_path);
    std::cout<<"basis set data read!"<<std::endl;
    total=0;
    for(int s=0;s!=n_sym;s++)
    {
       for(int t=0;t!=basis_size_sym[s];t++)
       {
          contraction_coeff[total]=new double[contraction_number_sym[s][t]];
          contraction_zeta[total]=new double[contraction_number_sym[s][t]];
          contraction_coeff_sym[s][t]=new double [contraction_number_sym[s][t]];
          contraction_zeta_sym[s][t]=new double [contraction_number_sym[s][t]];
          angular_mom_numbers_sym[s][t]=new int [2];
          angular_mom_numbers[total]=new int [2];
//          std::cout<<contraction_number_sym[s][t]<<" , "<<s<<" , "<<t<<std::endl;
          total++;
       }
    }
    std::cout<<"INITIALIZED ALL BASIS SET ARRAYS"<<std::endl;
    basis_data_reader(n_sym,basis_size_sym,contraction_number_sym,contraction_coeff_sym,contraction_zeta_sym,nucl_basis_func_sym,basis_func_type_sym,molpro_output_path,angular_mom_numbers_sym);
    total=0;
    int total2(0);
    for(int s=0;s!=n_sym;s++)
    {
       for(int t=0;t!=basis_size_sym[s];t++)
       {
          contraction_number[total]=contraction_number_sym[s][t];

          if(contraction_number[total] > max_contraction_num)
             max_contraction_num=contraction_number[total];

          basis_func_type[total]=basis_func_type_sym[s][t];
          nucl_basis_func[total]=nucl_basis_func_sym[s][t];
          for(int k=0;k!=contraction_number_sym[s][t];k++)
          {
             contraction_zeta[total][k]=contraction_zeta_sym[s][t][k];
             contraction_coeff[total][k]=contraction_coeff_sym[s][t][k];
          }
          angular_mom_numbers[total][0]=l_number(basis_func_type[total].c_str());
          angular_mom_numbers[total][1]=ml_number(basis_func_type[total].c_str(),angular_mom_numbers[total][0]);
//          std::cout<<basis_func_type[total]<<" !!"<<std::endl;
          total++;
       }
    }
    contraction_coeff_array=new double *[basis_size];
    contraction_zeta_array=new double *[basis_size];
    for(int i=0;i!=basis_size;i++)
    {
//       std::cout<<"nucl_basis_func "<<i<<" == "<<nucl_basis_func[i]<<std::endl;
       contraction_coeff_array[i]=new double[max_contraction_num];
       contraction_zeta_array[i]=new double[max_contraction_num];

       for(int j=0;j!=max_contraction_num;j++)
       {
          if(j<contraction_number[i])
          {
            contraction_coeff_array[i][j]=contraction_coeff[i][j];
            contraction_zeta_array[i][j]=contraction_zeta[i][j];
          }
          else
          {
            //contraction_coeff_array[i][j]=NAN;
            //contraction_zeta_array[i][j]=NAN;
            contraction_coeff_array[i][j]=0.0;
            contraction_zeta_array[i][j]=0.0;
          }
       }
    }
 }
    //*****************************THIS SECTION OF THE CODE DEPENDS ON NUCLEAR POSITIONS.*****************************
    std::cout<<"position "<<x<<std::endl;
    nucl_spher_pos[x][0][0]=rH/0.529177249;
    nucl_spher_pos[x][0][1]=acos(-1);
    nucl_spher_pos[x][0][2]=0.0;
    nucl_spher_pos[x][1][0]=rLi/0.529177249;
    nucl_spher_pos[x][1][1]=0.0;
    nucl_spher_pos[x][1][2]=0.0;
   std::cout<<"spher 1: "<<nucl_spher_pos[x][0][0]<<","<<nucl_spher_pos[x][0][1]<<","<<nucl_spher_pos[x][0][2]<<std::endl;
   std::cout<<"spher 2: "<<nucl_spher_pos[x][1][0]<<","<<nucl_spher_pos[x][1][1]<<","<<nucl_spher_pos[x][1][2]<<std::endl;
//    std::cout<<nucl_spher_pos[x][0][0]<<";"<<nucl_spher_pos[x][1][0]<<std::endl;

    //GET THE SIZE OF THE CI VECTOR IN THE NEUTRAL AND THE CATION
    if(symmetry)
    {
       num_of_ci_reader(n_states_neutral_sym, n_states_cat_sym, ci_size_neut_sym, ci_size_cat_sym, molpro_output_path,n_occs,n_sym);
       ci_size_neut=0;
       ci_size_cat=0;
       for(int i=0;i!=n_sym;i++)
       {
         ci_size_neut+=ci_size_neut_sym[i];
         std::cout<<"CI SIZE NEUT SYM "<<i+1<<" IS "<<ci_size_neut_sym[i]<<std::endl;
         ci_size_cat+=ci_size_cat_sym[i];
         std::cout<<"CI SIZE CAT SYM "<<i+1<<" IS "<<ci_size_cat_sym[i]<<std::endl;
       }
      // std::cout<<"PROBE CI SIZE QUERY DONE"<<std::endl;//DEBOGAGE
    }
    else
       num_of_ci_reader(&n_states_neut, &n_states_cat, &ci_size_neut, &ci_size_cat, molpro_output_path,n_occs);
    //std::cout<<"ci size neut is "<<ci_size_neut<<std::endl<<"ci size cat is "<<ci_size_cat<<std::endl;

    //ALLOCATE ARRAYS THAT DEPEND ON THE SIZE OF THE CI VECTOR
    ci_vec_neut[0]=new double[n_elec_neut*ci_size_neut+n_states_neut*ci_size_neut];//vector partitionned in two sections. Section 1 is filled with the mo label of each electron. section 2 is filled with CI coeff.
    ci_vec_neut[1]=new double[n_elec_neut*ci_size_neut];//this vector represents the spin state of each electron
    ci_vec_cat[0]=new double [(n_elec_neut-1)*ci_size_cat+n_states_cat*ci_size_cat];
    ci_vec_cat[1]=new double [(n_elec_neut-1)*ci_size_cat];

/*   for(int j=0;j!=num_of_nucl;j++)
   {
   std::cout<<"pos "<<x<<"spher "<<j<<": "<<nucl_spher_pos[x][j][0]<<","<<nucl_spher_pos[x][j][1]<<","<<nucl_spher_pos[x][j][2]<<std::endl;
   }
*/


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
//DEBOGAGE
/*
    for(int i=0;i!=ci_size_neut;i++)
    {
       for(int j=0;j!=n_elec_neut;j++)
       {
          std::cout<<"Neut conf "<<i<<" electron "<<j<<" => "<<int(ci_vec_neut[0][i*(n_elec_neut+n_states_neut)+j])<<std::endl;
       }
       for(int j=0;j!=n_states_neut;j++)
       {
          std::cout<<"Neut conf "<<i<<" state "<<j<<" => "<<double(ci_vec_neut[0][i*(n_elec_neut+n_states_neut)+j+n_elec_neut])<<std::endl;
       }
    }
    for(int i=0;i!=ci_size_cat;i++)
    {
       for(int j=0;j!=n_elec_neut-1;j++)
       {
          std::cout<<"cat conf "<<i<<" electron "<<j<<" => "<<int(ci_vec_cat[0][i*(n_elec_neut-1+n_states_cat)+j])<<std::endl;
       }
       for(int j=0;j!=n_states_cat;j++)
       {
          std::cout<<"cat conf "<<i<<" state "<<j<<" => "<<double(ci_vec_cat[0][i*(n_elec_neut-1+n_states_cat)+j+n_elec_neut-1])<<std::endl;
       }
    }
    exit(EXIT_SUCCESS);
    */
//DEBOGAGE
    MO_coeff_neutral[x]=new double[(n_occ+n_closed)*basis_size];
    overlap=new double[n_occ*n_occ];
    mo_dipole[x]=new double *[3];
    mo_dipole_spher[x]=new double *[3];
    mo_dipole[x][0]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole[x][1]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole[x][2]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[x][0]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[x][1]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[x][2]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    tran_den_mat_mo[x]=new double*[n_states_neut*n_states_neut];
    for(int i=0;i!=n_states_neut*n_states_neut;i++)
    {
      tran_den_mat_mo[x][i]=new double [(n_occ+n_closed)*(n_occ+n_closed) ];
    }

    std::cout<<"ALLOCATION OF CAS DEPENDENT VECTORS DONE"<<std::endl;
//*****************************COMPUTE DYSON ORBITALS*****************************
    //COMPUTE THE OVERLAP MATRIX BETWEEN THE MO OF THE NEUTRAL AND THE MO OF THE CATION FROM THE AO OVERLAP MATRIX
    if(symmetry)
       overlap_MO(overlap,n_occs,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral[x],n_sym);
    else
       overlap_MO(overlap,&n_occ,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral[x],n_sym);
    std::cout<<"MOLECULAR ORBITALS OVERLAP MATRIX READ WITHOUT ISSUE"<<std::endl;
    //COMPUTE THE MOLECULAR ORBITALS TRANSITION DIPOLE MOMENT MATRIX


    if(symmetry)
       dipole_MO(mo_dipole[x],n_occs,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    else
       dipole_MO(mo_dipole[x],&n_occ,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    std::cout<<"MOLECULAR ORBITALS DIPOLE MATRIX READ WITHOUT ISSUE"<<std::endl;

    for(int i=0;i!=(n_occ+n_closed)*(n_occ+n_closed);i++)
    {
       mo_dipole_spher[x][0][i]=sqrt(mo_dipole[x][0][i]*mo_dipole[x][0][i]+mo_dipole[x][1][i]*mo_dipole[x][1][i]+mo_dipole[x][2][i]*mo_dipole[x][2][i]);
       if(mo_dipole_spher[x][0][i]==0)
       {
          mo_dipole_spher[x][1][i]=0;
          mo_dipole_spher[x][2][i]=0;
       }
       else if(mo_dipole[x][0][i]==0)
       {
          mo_dipole_spher[x][1][i]=acos(mo_dipole[x][2][i]/mo_dipole_spher[x][0][i]);
          if(mo_dipole[x][1][i]==0)
          {
             mo_dipole_spher[x][2][i]=0;
          }
          else if (mo_dipole[x][1][i]>0)
          {
             mo_dipole_spher[x][2][i]=acos(-1)/2;
          }
          else
          {
             mo_dipole_spher[x][2][i]=3*acos(-1)/2;
          }
       }
       else
       {
           mo_dipole_spher[x][1][i]=acos(mo_dipole[x][2][i]/mo_dipole_spher[x][0][i]);
           mo_dipole_spher[x][2][i]=atan(mo_dipole[x][1][i]/mo_dipole[x][0][i]);
       }
//       std::cout<<"["<<int(i/n_occ)<<","<<int(i%n_occ)<<"]"<<"("<<mo_dipole[0][i]<<","<<mo_dipole[1][i]<<","<<mo_dipole[2][i]<<") => "<<"("<<mo_dipole_spher[0][i]<<","<<mo_dipole_spher[1][i]<<","<<mo_dipole_spher[2][i]<<")"<<std::endl;
    }
    //COMPUTE THE DYSON MO COEFFICIENTS IN THE BASIS OF THE MO OF THE NEUTRAL
    std::cout<<"COMPUTING DYSON MO COEFFICIENTS"<<std::endl;
//    DEBOGAGE!!!
    dyson_mo_basis_coeff[x]=new double[n_occ*n_states_neut*n_states_cat];
    dyson_mo_coeff_comp( n_states_neut,n_states_cat, n_occ,ci_size_neut, ci_size_cat, n_elec_neut, ci_vec_neut, ci_vec_cat,overlap, dyson_mo_basis_coeff[x]);

    std::cout<<"ENTERING DENSITY ROUTINE"<<std::endl;

    build_transition_density_matrix(n_states_neut,n_closed,n_occ,ci_size_neut,n_elec_neut,ci_vec_neut,tran_den_mat_mo[x]);

    //DEBOGAGE
   /* 
    double  valx(0);
    double  valy(0);
    double  valz(0);
    for(int inistate=0;inistate!=n_states_neut;inistate++)
    {
    for(int finstate=0;finstate!=n_states_neut;finstate++)
    {
    valx=0;
    valy=0;
    valz=0;
    for(int i = 0 ; i!=n_occ+n_closed ; i++)
    {
        for(int j = 0 ; j!=n_occ+n_closed ; j++)
        {
           valx+=tran_den_mat_mo[x][inistate*n_states_neut+finstate][i*(n_occ+n_closed)+j]*mo_dipole[x][0][i*(n_occ+n_closed)+j];
           valy+=tran_den_mat_mo[x][inistate*n_states_neut+finstate][i*(n_occ+n_closed)+j]*mo_dipole[x][1][i*(n_occ+n_closed)+j];
           valz+=tran_den_mat_mo[x][inistate*n_states_neut+finstate][i*(n_occ+n_closed)+j]*mo_dipole[x][2][i*(n_occ+n_closed)+j];
        }
    }
    if(inistate == finstate)
        valz+=(1.63/.529)*(3.*mH-1.*mLi)/(mLi+mH);
    std::cout<<"transition dipole "<<inistate<<"-"<<finstate<<":"<<valx<<","<<valy<<","<<valz<<std::endl;
    }
    }
//    exit(EXIT_SUCCESS);
//    */
    //DEBOGAGE

    std::cout<<" DENSITY ROUTINE DONE !"<<std::endl;

   }

for(int i=0;i!=grid_size;i++)
{
   for(int j=0;j!=num_of_nucl;j++)
   {
   std::cout<<"pos "<<i<<"spher "<<j<<": "<<nucl_spher_pos[i][j][0]<<","<<nucl_spher_pos[i][j][1]<<","<<nucl_spher_pos[i][j][2]<<std::endl;
   }
}
std::cout<<"POSITION DEPENDENT PART DONE"<<std::endl;
//    exit(EXIT_SUCCESS);
 
/*
 * BUILDING AND TESTING ONE ELECTRON REDUCED DENSITY MATRIX IN THE BASIS SET OF MO'S
 */

/*
for(int i = 0 ; i != (n_occ+n_closed)  ; i++)
{
   for(int j = 0 ; j != (n_occ+n_closed)  ; j++)
   {
      std::cout<<setw(10)<<setprecision(5)<<i<<","<<j<<" - "<<tran_den_mat_mo[0][0][i*(n_occ+n_closed)+j]<<std::endl;
   }
}

std::cout<<"********"<<std::endl;
for(int i = 0 ; i != (n_occ+n_closed)  ; i++)
{
   for(int j = 0 ; j != (n_occ+n_closed)  ; j++)
   {
      std::cout<<setw(10)<<setprecision(5)<<i<<","<<j<<" - "<<tran_den_mat_mo[0][1][i*(n_occ+n_closed)+j]<<std::endl;
   }
}
std::cout<<"********"<<std::endl;
for(int i = 0 ; i != (n_occ+n_closed)  ; i++)
{
   for(int j = 0 ; j != (n_occ+n_closed)  ; j++)
   {
      std::cout<<setw(10)<<setprecision(5)<<i<<","<<j<<" - "<<tran_den_mat_mo[0][2][i*(n_occ+n_closed)+j]<<std::endl;
   }
}
std::cout<<"********"<<std::endl;
*/
//exit(EXIT_SUCCESS);
/*
 * END OF DENSITY TESTING 
 */

//   build_ao_s(NULL,nucl_basis_func,contraction_number,nucl_spher_pos[0],contraction_coeff,contraction_zeta,basis_func_type,basis_size); 
/*
    double temp_norm;
    for(int i=0;i!=n_states_neut;i++)
    {
       for(int j=0;j!=n_states_cat;j++)
       {
          temp_norm=0;
          for(int t=0;t!=n_occ;t++)
          {
             temp_norm+=dyson_mo_basis_coeff[0][i*n_states_cat*n_occ+j*n_occ+t]*dyson_mo_basis_coeff[0][i*n_states_cat*n_occ+j*n_occ+t];
          }
          std::cout<<"Dyson orbital norm between states "<<i<<" and "<<j<<" : "<<temp_norm<<std::endl;
       }
    }
*/
/*
 * BUILDING AND PRINTING DYSON MO'S CUBES FOR TESTING 
 * */
/*
    std::string dyson_cube_loc("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/LiH_anal_dyson_");
    int nx(150);
    int ny(150);
    int nz(150);
    double cxmin(-37.5);
    double cxmax(-37.5+nx*0.5);
    double cymin(-37.5);
    double cymax(-37.5+nx*0.5);
    double czmin(-37.5);
    double czmax(-37.5+nx*0.5);
    int j(0);
    for(int i=0;i!=n_states_neut;i++)
    {
//       for(int j=0;j!=n_states_cat;j++)
       {
          std::cout<<"WRITING CUBE FOR DYSON "<<i<<" - "<<j<<std::endl;
         cube_header(&n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&num_of_nucl,&basis_size,contraction_number,nucl_coord[0],nucl_spher_pos[0],MO_coeff_neutral[0],dyson_mo_basis_coeff[0],contraction_coeff_array,contraction_zeta_array,nucl_basis_func,basis_func_type,angular_mom_numbers, dyson_cube_loc,i,j,2,nx, ny, nz, cxmin, cxmax, cymin, cymax, czmin, czmax);
       }
    }
      exit(EXIT_SUCCESS);
      */
/*
 * BUILDING AND PRINTING DYSON MO'S CUBES FOR TESTING 
 * */

    /*
     TESTING HF5 DIALOG
     */
      write_output(hf5_outfile, &n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&num_of_nucl,&basis_size,contraction_number,nucl_coord,nucl_spher_pos,mo_dipole,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_coeff_array,contraction_zeta_array,nucl_basis_func,basis_func_type,tran_den_mat_mo);
   //   read_output(hf5_outfile, &n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&num_of_nucl,&basis_size,nucl_spher_pos,mo_dipole,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_number,contraction_coeff,contraction_zeta,nucl_basis_func,basis_func_type);
    /*
     TESTING HF5 DIALOG
     */

    double kp(0);
    double kmax(1.5);
    int state_neut(0);
    int nk=256;
    int jl_max(6);
    int x(0);
    double *temp=new double[3];

    ifstream distrib_file;
    
    //We compute the values of the integral for every value of l1 and m1 ( = basis functions )  and for each value of l2 and m2 (Bessel function)
    //These integrals are also computed for the values of l3 going from 0 to l1+l2+1 and for m3 in -l3,l3
    //The integral is made of two factors: one depending on k, the other independent of k
    double begin=omp_get_wtime();

    //Then, we combine using LCAO coefficients to get MO's
    double*** pice_ortho_mo = new double **[jl_max*jl_max+2*jl_max+1];
    double*** pice_ddx_mo = new double **[jl_max*jl_max+2*jl_max+1];
    double*** pice_ddy_mo = new double **[jl_max*jl_max+2*jl_max+1];
    double*** pice_ddz_mo = new double **[jl_max*jl_max+2*jl_max+1];

    double ***pice_x = new double **[n_states_neut*n_states_cat];
    double ***pice_y = new double **[n_states_neut*n_states_cat];
    double ***pice_z = new double **[n_states_neut*n_states_cat];

   for(int mm=0;mm!=n_states_neut*n_states_cat;mm++)
   {
      pice_x[mm]=new double *[jl_max*jl_max+2*jl_max+1];
      pice_y[mm]=new double *[jl_max*jl_max+2*jl_max+1];
      pice_z[mm]=new double *[jl_max*jl_max+2*jl_max+1];
      for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
      {
         pice_x[mm][ji]=new double [nk];
         pice_y[mm][ji]=new double [nk];
         pice_z[mm][ji]=new double [nk];
      }
   }
   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   {
      pice_ortho_mo[ji] =  new double *[n_occ];
      pice_ddx_mo[ji] =  new double *[n_occ];
      pice_ddy_mo[ji] =  new double *[n_occ];
      pice_ddz_mo[ji] =  new double *[n_occ];

      for(int mm=0;mm!=n_occ;mm++)
      {
         pice_ortho_mo[ji][mm]=new double [nk];
         pice_ddx_mo[ji][mm]=new double [nk];
         pice_ddy_mo[ji][mm]=new double [nk];
         pice_ddz_mo[ji][mm]=new double [nk];
      }
   }
   compute_bessel_pice_mo(pice_ortho_mo,pice_ddx_mo,pice_ddy_mo,pice_ddz_mo,jl_max,n_occ,basis_size,nk,kmax,MO_coeff_neutral[x],contraction_zeta,contraction_coeff,contraction_number,nucl_spher_pos[x],nucl_basis_func,angular_mom_numbers);

   for(int nn=0;nn!=n_states_neut;nn++)
   {
      for( int nc=0;nc!=n_states_cat;nc++)
      {
         for(int mm=0;mm!=n_occ;mm++)
         {
            if(dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]!=0)
            {
               for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
               {
                  for(int k=0;k!=nk;k++)
                  {
                     pice_x[nn*n_states_cat+nc][ji][k]+=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*pice_ddx_mo[ji][mm][k];
                     pice_y[nn*n_states_cat+nc][ji][k]+=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*pice_ddy_mo[ji][mm][k];
                     pice_z[nn*n_states_cat+nc][ji][k]+=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*pice_ddz_mo[ji][mm][k];

                     for(int mmp=0;mmp!=n_occ;mmp++)
                     {
                        pice_x[nn*n_states_cat+nc][ji][k]-=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*mo_dipole[x][0][mm*n_occ+mmp]*pice_ortho_mo[ji][mmp][k];
                        pice_y[nn*n_states_cat+nc][ji][k]-=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*mo_dipole[x][1][mm*n_occ+mmp]*pice_ortho_mo[ji][mmp][k];
                        pice_z[nn*n_states_cat+nc][ji][k]-=dyson_mo_basis_coeff[x][nn*n_states_cat*n_occ+nc*n_occ+mm]*mo_dipole[x][2][mm*n_occ+mmp]*pice_ortho_mo[ji][mmp][k];
                        
                     }
                  }
               }
            }
         }
      }
   }
//   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   std::complex<double> ddtemp(0);
   for(int k=0;k!=nk;k++)
   {
      kp=(k+1)*kmax/nk;
      ddtemp=0;
//      std::cout<<kp<<","<<pice_z[0*n_states_cat][0][k]<<std::endl;
//      std::cout<<kp<<","<<pow(std::abs(pice_z[0*n_states_cat][0][k]),2)<<std::endl;
   for(int ji=0;ji!=jl_max*jl_max+2*jl_max+1;ji++)
   {
      ddtemp+=abs(pice_z[0*n_states_cat][ji][k]);
//      ddtemp+=(pice_z[0*n_states_cat][ji][k]);
   }
   std::cout<<kp*kp*27.211/2<<","<<kp*kp*pow(std::abs(ddtemp),2)<<std::endl;
   }


    double end=omp_get_wtime(); 
    std::cout<<std::endl<<std::fixed<<std::setprecision(6)<<(end-begin)/10<<"s is the time with a precision of"<<omp_get_wtick()<<std::endl;
       
//          std::cout<<kp*kp*27.211/2<<","<<int_cs<<std::endl;
//          
      return 0;
}

