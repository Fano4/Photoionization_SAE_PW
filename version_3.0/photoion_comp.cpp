#include "photoion_comp.hpp"

int main(int argc,char* argv[])
{
   int photoion_comp(int argc, char* argv[]);
   omp_set_num_threads(16); 

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
    const std::string hf5_outfile("LiH_1.6125_test.h5");
    int grid_size(1);
    double ***nucl_spher_pos=new double**[grid_size];
    double **MO_coeff_neutral=new double*[grid_size];
    double *overlap;
    double ***mo_dipole=new double**[grid_size];
    double ***mo_dipole_spher=new double**[grid_size];
    double **dyson_mo_basis_coeff=new double*[grid_size];
    double xp(0);
    double rLi(0);
    double rH(0);
    double *nucl_coord=new double[grid_size];

    //TEMPORARY VARIABLES
    double sum(0);
    int index(0);
    int index2(0);
    int temp_int(0);
    double xmin(1.6125);
    double xmax(1.6125);
    double mLi(6.015122795);
    double mH(1.007825);
    std::string file_root("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/LiH_");
    stringstream ss_molpro_file;
    std::string molpro_output_path;


//*****************************GET DATA FOR COMPUTING DYSON ORBITALS FROM MOLPRO OUTPUT FILE*****************************
        //GET THE NUMBER OF ELECTRONIC STATES AND THE SIZE OF THE ACTIVE SPACE
for(int x=0;x!=grid_size;x++)
{
    //*****************************THIS SECTION OF THE CODE DOES NOT DEPEND ON NUCLEAR POSITIONS. IT MUST BE CALLED ONLY ONCE*****************************
   //
   //

   xp=xmin+x*(xmax-xmin)/grid_size;
   rH=mLi*xp/(mLi+mH);
   rLi=mH*xp/(mLi+mH);

   nucl_coord[x]=xp/.529;
   ss_molpro_file.str("");
   ss_molpro_file<<file_root.c_str()<<xp<<".out";
   molpro_output_path=ss_molpro_file.str();

   std::cout<<"WORKING WITH MOLPRO FILE \""<<molpro_output_path.c_str()<<"\""<<std::endl;

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
    int max_contraction_num(0);
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
            contraction_coeff_array[i][j]=NAN;
            contraction_zeta_array[i][j]=NAN;
          }
       }
    }
   }
    //*****************************THIS SECTION OF THE CODE DEPENDS ON NUCLEAR POSITIONS.*****************************
    nucl_spher_pos[x]=new double*[num_of_nucl];
    for(int i=0;i!=num_of_nucl;i++)
    {
       nucl_spher_pos[x][i]=new double[3];
    }
    nucl_spher_pos[x][0][0]=rH/.529;
    nucl_spher_pos[x][0][1]=acos(-1);
    nucl_spher_pos[x][0][2]=0.0;
    nucl_spher_pos[x][1][0]=rLi/.529;
    nucl_spher_pos[x][1][1]=0.0;
    nucl_spher_pos[x][1][2]=0.0;

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

    dyson_mo_basis_coeff[x]=new double[n_occ*n_states_neut*n_states_cat];


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

    std::cout<<"ALLOCATION OF CAS DEPENDENT VECTORS DONE"<<std::endl;
//*****************************COMPUTE DYSON ORBITALS*****************************
    //COMPUTE THE OVERLAP MATRIX BETWEEN THE MO OF THE NEUTRAL AND THE MO OF THE CATION FROM THE AO OVERLAP MATRIX
    if(symmetry)
       overlap_MO(overlap,n_occs,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral[x],n_sym);
    else
       overlap_MO(overlap,&n_occ,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral[x],n_sym);
    std::cout<<"MOLECULAR ORBITALS OVERLAP MATRIX READ WITHOUT ISSUE"<<std::endl;
/*
    for(int i=0;i!=n_occ;i++)
    {
       for(int j=0;j!=basis_size;j++)
       {
          std::cout<<"geom "<<x<<" mo "<<i<<" basis function "<<j<<" , "<<MO_coeff_neutral[x][i*basis_size+j]<<std::endl;
       }
    }
    exit(EXIT_SUCCESS);
    */
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
    dyson_mo_coeff_comp( n_states_neut,n_states_cat, n_occ,ci_size_neut, ci_size_cat, n_elec_neut, ci_vec_neut, ci_vec_cat,overlap, dyson_mo_basis_coeff[x]);
/*    for(int i=0;i!=n_occ;i++)
    {
       for(int j=0;j!=basis_size;j++)
       {
          MO_coeff_neutral[x][i*basis_size+j]=0;
       }
          MO_coeff_neutral[x][i*basis_size+37]=1;
    }
    
      for(int i=0;i!=n_states_neut;i++)
      {
         for(int j=0;j!=n_occ;j++)
         {
            dyson_mo_basis_coeff[x][i*n_states_cat*n_occ+j]=0;
         }
         dyson_mo_basis_coeff[x][i*n_states_cat*n_occ+i]=1;
      }
//     DEBOGAGE!!!*/

    delete [] ci_vec_neut;
    delete [] ci_vec_cat;
    delete [] overlap;
   }

std::cout<<"POSITION DEPENDENT PART DONE"<<std::endl;
//    exit(EXIT_SUCCESS);
 

//   build_ao_s(NULL,nucl_basis_func,contraction_number,nucl_spher_pos[0],contraction_coeff,contraction_zeta,basis_func_type,basis_size); 

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
      //write_output(hf5_outfile, &n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&num_of_nucl,&basis_size,contraction_number,nucl_coord,nucl_spher_pos,mo_dipole,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_coeff_array,contraction_zeta_array,nucl_basis_func,basis_func_type);
//      read_output("test.h5", &n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&num_of_nucl,&basis_size,nucl_spher_pos,mo_dipole,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_number,contraction_coeff,contraction_zeta,nucl_basis_func,basis_func_type);
    /*
     TESTING HF5 DIALOG
     */

    double kp(0);
    double kmax(sqrt(2*200/27.211)/4);
    int state_neut(11);
    double int_cs(0);
    double thet(0.3);
    double phi(1.25);
    int nk=50;
    std::complex<double> *temp=new std::complex<double>[3];
//    std::complex<double> temp;
    std::complex<double> sinet;
    std::complex<double> cosinet;
    std::complex<double> sinep;
    std::complex<double> cosinep;
    std::complex<double> modulus;
    int n_theta=20;
    int n_phi=40;
    //double efield[3]={1,0,0};//SPHERICAL COORDINATES
    double efield[3]={0,0,0.05};//CARTESIAN COORDINATES
    int i(0);
    int t(0);
    int p(0);
    int j(0);
//    int basis_fun_index(7);
//    int contraction_index(1);


/*    temp=0;
    for(int k=0;k!=nk;k++)
{
   kp=k*kmax/nk;
       for( t=0;t!=n_theta;t++)
       {
          thet=t*(acos(-1))/n_theta;
          for( p=0;p!=n_phi;p++)
          {
             phi=2*p*acos(-1)/n_phi;

            temp+=pow(abs(MO_Fourier_transform(basis_fun_index,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)),2)*pow(kp,2)*sin(thet)*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)*(kmax/nk);
          }
       }
            std::cout<<temp<<std::endl;
}
std::cout<<"=>"<<real(temp)<<std::endl;
exit(EXIT_SUCCESS);
*/
    for(int k=0;k!=nk;k++)
    {
       kp=kmax*(k+1)/nk;
//       kp=0;
       int_cs=0;
//#pragma omp parallel for reduction(+:int_cs) private (i,j,t,p,temp,thet,phi)
       for( t=0;t<n_theta;t++)
       {
//          std::cout<<t<<std::endl;
          thet=t*(acos(-1))/n_theta;
          for( p=0;p<n_phi;p++)
          {
//            std::cout<<"point!"<<k<<","<<t<<","<<p<<"=>"<<temp<<std::endl;
             phi=2*p*acos(-1)/n_phi;
             temp[0]=0;
             temp[1]=0;
             temp[2]=0;
            // temp=0;
             for( i=0;i<n_occ;i++)
             {
                // df/dx = df/dk dk/dx + df/dthet dthet/dx + df/dphi dphi/dx 
                // df/dy = df/dk dk/dy + df/dthet dthet/dy + df/dphi dphi/dy 
                // df/dz = df/dk dk/dz + df/dthet dthet/dz + df/dphi dphi/dz 
                
                if(dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i] != 0)
                {
                   temp[0]-=std::complex<double>(0,1)*dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *(
                         sin(thet) * cos(phi) * MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                         + cos(thet) * cos(phi) * MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                         - sin(phi) * MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                      );
/*                   std::cout<<kp<<","<<thet<<","<<phi<<" ; "<<i<<":"<<std::endl<<
                   "dyson "<<dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i]<<std::endl<<
                   "1: "<<std::complex<double>(0,1)*sin(thet) * cos(phi) * MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)<<std::endl<<
                   "2: "<<std::complex<double>(0,1)*cos(thet) * cos(phi) * MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)<<std::endl<<
                   "3: "<<- std::complex<double>(0,1)*sin(phi) * MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)<<std::endl<<std::endl;
  */              
                   temp[1]-=std::complex<double>(0,1)*dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *(
                         sin(thet) * sin(phi) * MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                       + cos(thet) * sin(phi) * MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                       + cos(phi) * MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                       );

                   temp[2]-=std::complex<double>(0,1)*dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *(
                         cos(thet) * MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                       - sin(thet) * MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                       );
                }
   /*              
                 temp-=std::complex<double>(0,1)*dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+i]
                          *(
                                efield[0]*(
                  sin(thet)*cos(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                  +cos(thet)*cos(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
               -sin(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                                   )
                                +efield[1]*(
                  sin(thet)*sin(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                  +cos(thet)*sin(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                  +cos(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                                   )
                                +efield[2]*(
                  cos(thet)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                  -sin(thet)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)
                                   ));
          
                 
                for(int j=0;j!=n_occ;j++)
                {
                
                   temp-=dyson_mo_basis_coeff[0][state_neut*n_states_cat*n_occ+0*n_occ+j]
                      *MO_Fourier_transform(i,kp,thet,phi,nucl_spher_pos[0],nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,angular_mom_numbers,MO_coeff_neutral[0],basis_size)*(mo_dipole[0][0][i*n_occ+j]*efield[0]+mo_dipole[0][1][i*n_occ+j]*efield[1]+mo_dipole[0][2][i*n_occ+j]*efield[2]);
 //                  exit(EXIT_SUCCESS);
                 }
 //            */
             }
             std::cout<<kp<<"    "<<thet<<"    "<<phi<<"    "<<real(temp[0])<<"   "<<imag(temp[0])<<"  "<<real(temp[1])<<"  "<<imag(temp[1])<<"   "<<real(temp[2])<<"  "<<imag(temp[2])<<std::endl;
//             std::cout<<thet<<","<<phi<<","<<pow(abs(temp),2)*pow(kp,2)*sin(thet)*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)<<std::endl;
             

            // int_cs+=pow(kp,1)*sin(thet)*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)*abs(temp)*abs(temp);
         //std::cout<<kp<<" , "<<thet<<" , "<<phi<<" , "<<temp<<std::endl;
//             std::cout<<"small loop "<<p<<std::endl;
          }std::cout<<std::endl;
//             std::cout<<"LARGE loop "<<t<<std::endl<<"##################################################################"<<std::endl;
       }
       
//          std::cout<<kp*kp*27.211/2<<","<<int_cs<<std::endl;
          
    }

      return 0;
}

