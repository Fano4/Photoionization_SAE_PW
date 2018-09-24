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
    int n_states_neutral_sym[n_sym]={10,4,4,1};//{8,3,3,1};
    int n_states_cat(0);
    int n_states_cat_sym[n_sym]={2,1,1,0};//{2,1,1,0}
    int n_occ(0);
    int *n_occs;
    int n_closed(0);
    int *n_closeds;
    string mfpad_name("mfpad_CH4_0_0.txt");
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
    int ci_size_neut(0);
    int ci_size_cat(0);
    int ci_size_neut_sym[n_sym];
    int ci_size_cat_sym[n_sym];
    int n_elec_neut(4);//!!!! ecrire une routine qui cherche le nombre d'electrons dans l'output molpro!!!
   double Pi=acos(-1);


    //TEMPORARY VARIABLES
    double kmin(0.0);
    double Emax(2.756238286);
    double kmax(sqrt(2*Emax));
    int nk(150);
    int ntheta(90);
    int nphi(120);
    int nx(300);//125);
    int ny(300);//125);
    int nz(300);//125);
    int nkx(64);
    int nky(64);
    int nkz(64);
    double xmin(-75);//-27.401029);
    double xmax(-75.0+nx*0.5);//27.842971);
    double ymin(-75.0);//-27.401029);
    double ymax(-75.0+ny*0.5);//27.842971);
    double zmin(-75.0);//-27.010679);
    double zmax(-75.0+nz*0.5);//28.233321);
    double x;
    double y;
    double z;
    int n_points_sphere(128);
    int continuum_matrix_cols(128);
    double Redip_mom[3];
    double Imdip_mom[3];
    double *Reinput;
    double *Iminput;
    double **neut_mo_cube_array;
    double *mo_cube_array;
    double **moment_cube_array;
    double **coord;
    double *overlap;
    double *dyson_mo_basis_coeff;
    double *ci_vec_neut[2];
    double *ci_vec_cat[2];
    double **sphere_dist=new double*[2];
    

    double energy(0);
    double sum(0);
    int index(0);
    int index2(0);
//    double kp(0);
    double *theta=new double[continuum_matrix_cols];
    double *phi=new double[continuum_matrix_cols];
    double random(0);
    double random2(0);
    double random3(0);
    double temp(0);
    int temp_int(0);

    string MO_cube_loc("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/lih_neut_orbital_");//"MO_CUBE_LOC");
    string dyson_cube_loc("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/dyson_mo_");//"DYSON_CUBE_LOC");
    stringstream ss_cross_section;
    string s_cross_section;
    stringstream ss_PICE;
    string s_PICE;

//    kp=kmin+18*(kmax-kmin)/nk;
    srand(time(NULL));

    mo_cube_array=new double[nx*ny*nz];
    moment_cube_array=new double*[3];
    coord=new double*[3];
    coord[0]=new double[nx*ny*nz];
    coord[1]=new double[nx*ny*nz];
    coord[2]=new double[nx*ny*nz];

        for(int i=0;i!=nx;i++)
        {
           x=xmin+i*(xmax-xmin)/nx;
           for(int j=0;j!=ny;j++)
           {
             y=ymin+j*(ymax-ymin)/ny;
              for(int m=0;m!=nz;m++)
              {
                 z=zmin+m*(zmax-zmin)/nz;
                 coord[0][i*ny*nz+j*nz+m]=x;
                 coord[1][i*ny*nz+j*nz+m]=y;
                 coord[2][i*ny*nz+j*nz+m]=z;
              }
           }
        }
double *rectemp=new double [nx*ny*nz];
double *imctemp=new double [nx*ny*nz];
double *xim_in=new double[nx*ny*nz];
double *temp_cub=new double[nx*ny*nz];
MKL_LONG status;
double param(1/sqrt(nx*ny*nz));
MKL_LONG ft_dim[3];
ft_dim[0]=nx;
ft_dim[1]=ny;
ft_dim[2]=nz;
DFTI_DESCRIPTOR_HANDLE cube_ft_desc;
       std::cout<<"PROBE ARRAYS ALLOCATION DONE"<<std::endl;//DEBOGAGE
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
    //ALLOCATE ARRAYS THAT DEPEND ON n_occ
   overlap=new double[n_occ*n_occ];
   neut_mo_cube_array=new double*[n_occ];
   double **re_ft_mo_array=new double *[n_occ];
   double **im_ft_mo_array=new double *[n_occ];
   double ***moment_orbital=new double **[3];
   moment_orbital[0]=new double *[n_occ];
   moment_orbital[1]=new double *[n_occ];
   moment_orbital[2]=new double *[n_occ];
   double **gradx_re_ft_mo_array=new double *[n_occ];
   double **gradx_im_ft_mo_array=new double *[n_occ];
   double **grady_re_ft_mo_array=new double *[n_occ];
   double **grady_im_ft_mo_array=new double *[n_occ];
   double **gradz_re_ft_mo_array=new double *[n_occ];
   double **gradz_im_ft_mo_array=new double *[n_occ];
   double **transition_dipole=new double*[3];
   transition_dipole[0]=new double[n_occ*n_occ];
   transition_dipole[1]=new double[n_occ*n_occ];
   transition_dipole[2]=new double[n_occ*n_occ];

    for(int k=0;k!=n_occ;k++)
    {
       neut_mo_cube_array[k]=new double[nx*ny*nz];
    }
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
       /*for(int i=0;i!=ci_size_neut;i++)
       {
          for(int j=0;j!=n_elec_neut;j++)
             std::cout<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*i+j]<<" ";
          std::cout<<"  ";

          for(int j=0;j!=n_states_neut;j++)
             std::cout<<ci_vec_neut[0][(n_elec_neut+n_states_neut)*i+n_elec_neut+j]<<" ";
          std::cout<<std::endl;
          for(int j=0;j!=n_elec_neut;j++)
             std::cout<<ci_vec_neut[1][(n_elec_neut)*i+j]<<" ";
       std::cout<<endl;
     }//DEBOGAGE
       std::cout<<std::endl<<"CATION"<<std::endl;
       for(int i=0;i!=ci_size_cat;i++)
       {
          for(int j=0;j!=n_elec_neut-1;j++)
             std::cout<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*i+j]<<" ";
          std::cout<<"  ";

          for(int j=0;j!=n_states_cat;j++)
             std::cout<<ci_vec_cat[0][(n_elec_neut-1+n_states_cat)*i+n_elec_neut-1+j]<<" ";
          std::cout<<std::endl;
          for(int j=0;j!=n_elec_neut-1;j++)
             std::cout<<ci_vec_cat[1][(n_elec_neut-1)*i+j]<<" ";
       std::cout<<endl;
     }//DEBOGAGE*/
    std::cout<<"CI_VEC_READER ROUTINE ENDED WITHOUT ISSUE"<<std::endl;

   //GET THE CUBES OF THE MO OF THE NEUTRAL

   if(symmetry)
   {
   temp_int=0;
   for(int l=0;l!=n_sym;l++)
   {
      for(int k=0;k!=n_occs[l];k++)
      {
         std::cout<<"probe cube "<<k+1<<"."<<l+1<<" => "<<temp_int+1<<std::endl;//DEBOGAGE 
         cube_reader(k,l,nx,ny,nz,MO_cube_loc,neut_mo_cube_array[temp_int]);
         temp_int++;
      }
   }
   }
   else
   {
      for(int k=0;k!=n_occ;k++)
      {
         std::cout<<"probe cube "<<k+1<<"."<<1<<" => "<<temp_int+1<<std::endl;//DEBOGAGE 
         cube_reader(k,0,nx,ny,nz,MO_cube_loc,neut_mo_cube_array[temp_int]);
         temp_int++;
      }
   }
    std::cout<<"MO CUBE READER ROUTINE ENDED WITHOUT ISSUE"<<std::endl;
//*****************************COMPUTE DYSON ORBITALS*****************************
    //COMPUTE THE OVERLAP MATRIX BETWEEN THE MO OF THE NEUTRAL AND THE MO OF THE CATION FROM THE AO OVERLAP MATRIX
    if(symmetry)
    overlap_MO(overlap,n_occs,&basis_size,molpro_output_path,n_sym);
    else
    overlap_MO(overlap,&n_occ,&basis_size,molpro_output_path,n_sym);
    /*for(int i=0;i!=n_occ;i++)
    {
       for(int j=0;j!=n_occ;j++)
       {
          std::cout<<setw(13)<<overlap[i*n_occ+j];
       }std::cout<<std::endl<<std::endl;
    }*/
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
       cube_header(dyson_mo_basis_coeff,n_occ,n_states_neut,n_states_cat,neut_mo_cube_array,dyson_cube_loc.c_str(),a,0,2,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,mo_cube_array);
    }

    exit(EXIT_SUCCESS);
    //BUILD THE CUBE OF THE DYSON ORBITALS
    //
std::cout<<"BEGINNING COMPUTATION OF PICE"<<std::endl;

status=DftiCreateDescriptor(&cube_ft_desc, DFTI_DOUBLE, DFTI_COMPLEX ,3, ft_dim);
std::cout<<"DFT initialization : "<<DftiErrorMessage(status)<<std::endl;
DftiSetValue(cube_ft_desc,DFTI_FORWARD_SCALE,param);
DftiSetValue(cube_ft_desc,DFTI_COMPLEX_STORAGE,DFTI_REAL_REAL);
DftiSetValue(cube_ft_desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
status=DftiCommitDescriptor(cube_ft_desc);
std::cout<<"Commiting DFT : "<<DftiErrorMessage(status)<<std::endl;
for(int k=0;k!=n_occ;k++)
{
   moment_orbital[0][k]=new double [nx*ny*nz];
   moment_orbital[1][k]=new double [nx*ny*nz];
   moment_orbital[2][k]=new double [nx*ny*nz];
   vdMul(nx*ny*nz,coord[0],neut_mo_cube_array[k],moment_orbital[0][k]);
   vdMul(nx*ny*nz,coord[1],neut_mo_cube_array[k],moment_orbital[1][k]);
   vdMul(nx*ny*nz,coord[2],neut_mo_cube_array[k],moment_orbital[2][k]);
}
for(int k=0;k!=n_occ;k++)
{
   for(int kp=0;kp!=n_occ;kp++)
   {
      transition_dipole[0][k*n_occ+kp]=cblas_ddot(nx*ny*nz,neut_mo_cube_array[k],1,moment_orbital[0][kp],1);
      transition_dipole[1][k*n_occ+kp]=cblas_ddot(nx*ny*nz,neut_mo_cube_array[k],1,moment_orbital[1][kp],1);
      transition_dipole[2][k*n_occ+kp]=cblas_ddot(nx*ny*nz,neut_mo_cube_array[k],1,moment_orbital[2][kp],1);
      //std::cout<<transition_dipole[2][k*n_occ+kp]<<"    ";
   }//std::cout<<std::endl;
}
for(int k=0;k!=n_occ;k++)
{
   std::cout<<"FOURIER TRANSFORM OF ORBITAL "<<k<<std::endl;
   im_ft_mo_array[k]=new double [nx*ny*nz];
   re_ft_mo_array[k]=new double [nx*ny*nz];
   gradx_re_ft_mo_array[k]=new double[nx*ny*nz];
   gradx_im_ft_mo_array[k]=new double[nx*ny*nz];
   grady_re_ft_mo_array[k]=new double[nx*ny*nz];
   grady_im_ft_mo_array[k]=new double[nx*ny*nz];
   gradz_re_ft_mo_array[k]=new double[nx*ny*nz];
   gradz_im_ft_mo_array[k]=new double[nx*ny*nz];
//STEP 1: FOURIER TRANSFORM THE MO'S
   /*for(int i=0;i!=nx*ny*nz;i++)
   {
      std::cout<<neut_mo_cube_array[k][i]<<std::endl;
   }*/
   cblas_dcopy(nx*ny*nz,neut_mo_cube_array[k],1,temp_cub,1);
   center_wave(temp_cub,ft_dim,3);
   cblas_dscal (nx*ny*nz, 0, xim_in, 1);
   std::cout<<"FT of MO "<<k<<": "<<DftiComputeForward(cube_ft_desc, temp_cub, xim_in, re_ft_mo_array[k],im_ft_mo_array[k])<<std::endl;
//   shift_phase(re_ft_mo_array[k],im_ft_mo_array[k],nx, ny, nz);
   center_wave(re_ft_mo_array[k],ft_dim,3);
   center_wave(im_ft_mo_array[k],ft_dim,3);
//STEP 2: COMPUTE THE GRADIENT OF THE MO'S IN THE RECIPROCAL SPACE
   center_wave(moment_orbital[0][k],ft_dim,3);
   cblas_dscal (nx*ny*nz, 0, xim_in, 1);
   std::cout<<" FT gradient along x of MO "<<k<<": "<<DftiComputeForward(cube_ft_desc,moment_orbital[0][k],xim_in,gradx_im_ft_mo_array[k],gradx_re_ft_mo_array[k])<<std::endl;
   cblas_dscal (nx*ny*nz, -1, gradx_im_ft_mo_array[k], 1);
   center_wave(gradx_re_ft_mo_array[k],ft_dim,3);
   center_wave(gradx_im_ft_mo_array[k],ft_dim,3);
   center_wave(moment_orbital[1][k],ft_dim,3);
   cblas_dscal (nx*ny*nz, 0, xim_in, 1);
   std::cout<<" FT gradient along y of MO "<<k<<": "<<DftiComputeForward(cube_ft_desc,moment_orbital[1][k],xim_in,grady_im_ft_mo_array[k],grady_re_ft_mo_array[k])<<std::endl;
   cblas_dscal (nx*ny*nz, -1, grady_im_ft_mo_array[k], 1);
   center_wave(grady_re_ft_mo_array[k],ft_dim,3);
   center_wave(grady_im_ft_mo_array[k],ft_dim,3);
   center_wave(moment_orbital[2][k],ft_dim,3);
   cblas_dscal (nx*ny*nz, 0, xim_in, 1);
   std::cout<<" FT gradient along z of MO "<<k<<": "<<DftiComputeForward(cube_ft_desc,moment_orbital[2][k],xim_in,gradz_im_ft_mo_array[k],gradz_re_ft_mo_array[k])<<std::endl;
   cblas_dscal (nx*ny*nz, -1, gradz_im_ft_mo_array[k], 1);
   center_wave(gradz_re_ft_mo_array[k],ft_dim,3);
   center_wave(gradz_im_ft_mo_array[k],ft_dim,3);
}
    
//*****************************COMPUTE IONIZATION MATRIX ELEMENTS IN THE ORTHOGONALIZED PLANE WAVE APPROX******************************USING THE PROPERTIES OF THE FOURIER TRANSFORM OF MOLECULAR ORBITALS

       ofstream PICE;
for(int m=0;m!=n_states_neut;m++)
{
   for(int n=0;n!=n_states_cat;n++)
   {
      std::cout<<"PRINTING PICE BETWEEN STATES "<<m<<","<<n<<std::endl;
      // #############################   X COMPONENT OF THE PICE VECTOR FIELD ##############
      cblas_dscal (nx*ny*nz, 0, rectemp, 1);
      cblas_dscal (nx*ny*nz, 0, imctemp, 1);
      for(int k=0;k!=n_occ;k++)
      {
         cblas_daxpy (nx*ny*nz,dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],gradx_im_ft_mo_array[k],1,rectemp,1);
         cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],gradx_re_ft_mo_array[k],1,imctemp,1);
         for(int kp=0;kp!=n_occ;kp++)
         {
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[0][k*n_occ+kp],re_ft_mo_array[kp],1,rectemp,1);
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[0][k*n_occ+kp],im_ft_mo_array[kp],1,imctemp,1);
         }
      }
     // cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),rectemp,1);
     // cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),imctemp,1);
     // #############################   REAL PART OF THE X COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_RePICE_X_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Real part, X component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nkx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<nky<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nkz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<rectemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
       PICE.close();
      // #############################   IMAGINARY PART OF THE X COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_ImPICE_X_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Imaginary part, X component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nkx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<nky<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nkz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<imctemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
       PICE.close();
      // #############################   Y COMPONENT OF THE PICE VECTOR FIELD ##############
     cblas_dscal (nx*ny*nz, 0, rectemp, 1);
     cblas_dscal (nx*ny*nz, 0, imctemp, 1);
      for(int k=0;k!=n_occ;k++)
      {
         cblas_daxpy (nx*ny*nz,dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],grady_im_ft_mo_array[k],1,rectemp,1);
         cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],grady_re_ft_mo_array[k],1,imctemp,1);
         for(int kp=0;kp!=n_occ;kp++)
         {
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[1][k*n_occ+kp],re_ft_mo_array[kp],1,rectemp,1);
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[1][k*n_occ+kp],im_ft_mo_array[kp],1,imctemp,1);
         }
      }
    //  cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),rectemp,1);
    //  cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),imctemp,1);
      // #############################   REAL PART OF THE Y COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_RePICE_Y_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Real part, Y component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<ny<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<rectemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
       PICE.close();
      // #############################   IMAGINARY PART OF THE Y COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_ImPICE_Y_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Imaginary part, Y component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nkx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<nky<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nkz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<imctemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
       PICE.close();
      // #############################   Z COMPONENT OF THE PICE VECTOR FIELD ##############
      cblas_dscal (nx*ny*nz, 0, rectemp, 1);
      cblas_dscal (nx*ny*nz, 0, imctemp, 1);
      for(int k=0;k!=n_occ;k++)
      {
         cblas_daxpy (nx*ny*nz,dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],gradz_im_ft_mo_array[k],1,rectemp,1);
         cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k],gradz_re_ft_mo_array[k],1,imctemp,1);
         for(int kp=0;kp!=n_occ;kp++)
         {
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[2][k*n_occ+kp],re_ft_mo_array[kp],1,rectemp,1);
            cblas_daxpy (nx*ny*nz,-dyson_mo_basis_coeff[m*n_occ*n_states_cat+n*n_occ+k]*transition_dipole[2][k*n_occ+kp],im_ft_mo_array[kp],1,imctemp,1);
         }
      }
   //   cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),rectemp,1);
   //   cblas_dscal(nx*ny*nz,sqrt(pow(2*acos(-1),3)/((xmax-xmin)*(ymax-ymin)*(zmax-zmin))),imctemp,1);
      // #############################   REAL PART OF THE Z COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_RePICE_Z_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Real part, Z component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nkx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<nky<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nkz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
       //integrated_cs=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<rectemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
        //std::cout<<"Total integrated cross section is "<<integrated_cs<<std::endl;
       PICE.close();
      // #############################   IMAGINARY PART OF THE Z COMPONENT ##############
       ss_PICE.str("");
       ss_PICE<<"LiH_ImPICE_Z_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();
       PICE.open(s_PICE.c_str());
       PICE<<"Photoionization coupling elements cube file"<<std::endl<<"Coupling between neutral state "<<m<<" and cation state "<<n<<"(Imaginary part, Z component ) \n";
       PICE<<1<<"  "<<-(acos(-1)*nkx/(xmax-xmin))<<"  "<<-(acos(-1)*nky/(ymax-ymin))<<"  "<<-(acos(-1)*nkz/(zmax-zmin))<<" \n";
       PICE<<nkx<<"   "<<2*acos(-1)/(xmax-xmin)<<"    0.000000    0.000000 \n";
       PICE<<nky<<"   0.000000    "<<2*acos(-1)/(ymax-ymin)<<"    0.000000 \n";
       PICE<<nkz<<"   0.000000    0.000000   "<<2*acos(-1)/(zmax-zmin)<<" \n";
       PICE<<"    1    1.000000    0.0000000000        0.0000000000      0.000000000 \n";
       PICE<<"1 111 \n";
       index=0;
       index2=0;
        for(int i=nx/2-nkx/2;i!=nx/2+nkx/2;i++)
        {
           for(int j=ny/2-nky/2;j!=ny/2+nky/2;j++)
           {
              for(int k=nz/2-nkz/2;k!=nz/2+nkz/2;k++)
              {
                PICE<<scientific<<setw(16)<<imctemp[i*ny*nz+j*nz+k];
                index++;
                index2++;
                if (index2%6==0)
                {
                    PICE<<"\n";
                }
              }
           }
        }
       PICE.close();
   }
}
//*****************************MAP THE COUPLING ELEMENTS OF THE CUBE ON A SPHERE IN THE RECIPROCAL SPACE EITHER WITH A RANDOM DISTRIBUTION OR WITH AN ORDERED GRID******************************
//
/*
for(int m=0;m!=n_states_neut;m++)
{
   for(int n=0;n!=n_states_cat;n++)
   {
       ss_PICE.str("");
       ss_PICE<<"LiH_PICE_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();

       PICE.open(s_PICE.c_str());

       for(int e=1;e!=nk+1;e++)
       {
          kp=e*kmax/nk;
          sum=0;
//   std::cout<<"The factor is "<<4*0.441*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)*kp/(2*Pi*137)<<std::endl;

          for(int i=0;i!=n_points_sphere;i++)
          {
            // std::cout<<"Entering loop"<<std::endl;
                theta[j]=sphere_dist[0][i*continuum_matrix_cols+j];
                phi[j]=sphere_dist[1][i*continuum_matrix_cols+j];
           //  std::cout<<"building pw..."<<std::endl;

             PICE<<kp<<"    "<<theta[w]<<"    "<<phi[w]<<"    "<<sqrt(2)*Reproduct_vector[0][w]<<"   "<<sqrt(2)*Improduct_vector[0][w]<<"  "<<sqrt(2)*Reproduct_vector[1][w]<<"  "<<sqrt(2)*Improduct_vector[1][w]<<"   "<<sqrt(2)*Reproduct_vector[2][w]<<"  "<<sqrt(2)*Improduct_vector[2][w]<<std::endl;
             if((i*continuum_matrix_cols+w+1)%nphi==0)
             {
                PICE<<std::endl;
             }
          }
       }
       PICE.close();
   }
}
*/
    return 0;
    
}
bool center_wave(double *x,MKL_LONG *grid_size, int num_of_dim)
{
   using namespace std;
//   double* temp=new double [grid_size[0]*grid_size[1]*grid_size[2]];

   for(int i=0;i!=grid_size[0];i++)
   {
      for(int j=0;j!=grid_size[1];j++)
      {
         for(int k=0;k!=grid_size[2]/2;k++)
         {
            swap(x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2],x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]);
//            x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2]=0.5*(x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2+1]+x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2-1]);
      //      temp[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]=x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2];
      //      temp[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k+grid_size[2]/2]=x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k];
         }
      }
   }
   for(int i=0;i!=grid_size[0];i++)
   {
      for(int j=0;j!=grid_size[1]/2;j++)
      {
         for(int k=0;k!=grid_size[2];k++)
         {
            swap(x[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2)*grid_size[2]+k],x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]);
       //     x[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2)*grid_size[2]+k]=0.5*(x[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2+1)*grid_size[2]+k]+x[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2-1)*grid_size[2]+k]);
     //       temp[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]=x[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2)*grid_size[2]+k];
     //       temp[i*grid_size[1]*grid_size[2]+(j+grid_size[1]/2)*grid_size[2]+k]=x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k];
         }
      }
   }
   for(int i=0;i!=grid_size[0]/2;i++)
   {
      for(int j=0;j!=grid_size[1];j++)
      {
         for(int k=0;k!=grid_size[2];k++)
         {
            swap(x[(i+grid_size[0]/2)*grid_size[1]*grid_size[2]+j*grid_size[2]+k],x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]);
    //        x[(i+grid_size[0]/2)*grid_size[1]*grid_size[2]+j*grid_size[2]+k]=0.5*(x[(i+grid_size[0]/2+1)*grid_size[1]*grid_size[2]+j*grid_size[2]+k]+x[(i+grid_size[0]/2-1)*grid_size[1]*grid_size[2]+j*grid_size[2]+k]);
            //temp[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k]=x[(i+grid_size[0]/2)*grid_size[1]*grid_size[2]+j*grid_size[2]+k];
            //temp[(i+grid_size[0]/2)*grid_size[1]*grid_size[2]+j*grid_size[2]+k]=x[i*grid_size[1]*grid_size[2]+j*grid_size[2]+k];
         }
      }
   }
 //  cblas_dcopy(grid_size[0]*grid_size[1]*grid_size[2],temp,1,x,1);
 //  delete [] temp;
    return 0;
}
/*void spherical_extract_from_cube(double k,double** sphere_dist,double n_points_sphere,double xmin,double xmax,double nx,double ymin,double ymax,double ny,double zmin,double zmax,double nz,double* Recube,double* Imcube,double* Resphere,double* Imsphere )
{
   double xsphere(0);
   double ysphere(0);
   double zsphere(0);
   int x_index(0);
   int y_index(0);
   int z_index(0);

   for(int i=0;i!=n_points_sphere;i++)
   {
      xsphere=k*sin(sphere_dist[0][i])*cos(sphere_dist[1][i]);
      ysphere=k*sin(sphere_dist[0][i])*sin(sphere_dist[1][i]);
      zsphere=k*cos(sphere_dist[0][i]);

      x_index=int((xsphere-xmin)*nx/(xmax-xmin));
      y_index=int((ysphere-ymin)*ny/(ymax-ymin));
      z_index=int((zsphere-zmin)*nz/(zmax-zmin));

      Resphere[i]=Recube[x_index*ny*nz+y_index*nz+z_index];
      Imsphere[i]=Imcube[x_index*ny*nz+y_index*nz+z_index];
   }
}*/
