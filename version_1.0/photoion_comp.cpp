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
    int n_states_neutral_sym[n_sym]={10,4,4,1};
    int n_states_cat(0);
    int n_states_cat_sym[n_sym]={4,1,1,0};//{2,1,1,0}
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
    double Emax(200/27.211);
    double kmax(sqrt(2*Emax)/4);
    int nk(50);
    int ntheta(20);
    int nphi(40);
    double xmin(-62.5);
    double xmax(62.5);
    double ymin(-62.5);
    double ymax(62.5);
    double zmin(-62.5);
    double zmax(62.5);
    int nx(180);
    int ny(180);
    int nz(180);
    double x;
    double y;
    double z;
    int n_points_sphere(ntheta*nphi);
    int continuum_matrix_cols(1);
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
    double kp(0);
    double *theta=new double[continuum_matrix_cols];
    double *phi=new double[continuum_matrix_cols];
    double random(0);
    double random2(0);
    double random3(0);
    double temp(0);
    int temp_int(0);
    double **Reproduct_vector=new double*[3];
    double **Improduct_vector=new double*[3];
    Reproduct_vector[0]=new double[continuum_matrix_cols];
    Reproduct_vector[1]=new double[continuum_matrix_cols];
    Reproduct_vector[2]=new double[continuum_matrix_cols];
    Improduct_vector[0]=new double[continuum_matrix_cols];
    Improduct_vector[1]=new double[continuum_matrix_cols];
    Improduct_vector[2]=new double[continuum_matrix_cols];

    string MO_cube_loc("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/lih_neut_orbital_");
    string dyson_cube_loc("/data1/home/stephan/LiH_gridtest_+++custom_MO_1.6125/dyson_mo_bruteforce");
    stringstream ss_cross_section;
    string s_cross_section;
    stringstream ss_PICE;
    string s_PICE;

    kp=kmin+18*(kmax-kmin)/nk;
    srand(time(NULL));

    Reinput=new double[continuum_matrix_cols*nx*ny*nz];
    Iminput=new double[continuum_matrix_cols*nx*ny*nz];
    mo_cube_array=new double[nx*ny*nz];
    moment_cube_array=new double*[3];
    moment_cube_array[0]=new double[nx*ny*nz];
    moment_cube_array[1]=new double[nx*ny*nz];
    moment_cube_array[2]=new double[nx*ny*nz];
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
       std::cout<<"PROBE ARRAYS ALLOCATION DONE"<<std::endl;//DEBOGAGE
//*****************************GET DATA FOR COMPUTING DYSON ORBITALS FROM MOLPRO OUTPUT FILE*****************************
        //GET THE NUMBER OF ELECTRONIC STATES AND THE SIZE OF THE ACTIVE SPACE
//    n_states_reader(&n_states_neut,&n_states_cat,&n_elec_neut,molpro_output_path);

       std::cout<<molpro_output_path<<std::endl;
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

    //BUILD THE CUBE OF THE DYSON ORBITALS
/*    for(int i=0;i!=n_states_neut;i++)
    {
        cube_header(dyson_mo_basis_coeff,n_occ,n_states_neut,n_states_cat,neut_mo_cube_array,dyson_cube_loc,i,0,num_of_nucl,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,mo_cube_array); 
    }
      std::cout<<"DYSON ORBITAL CUBE CONSTRUCTED"<<std::endl;
*/

//*****************************COMPUTE IONIZATION MATRIX ELEMENTS IN THE ORTHOGONALIZED PLANE WAVE APPROX******************************
    sphere_dist[0]=new double[n_points_sphere];
    sphere_dist[1]=new double[n_points_sphere];
    double thet(0);
    double php(0);
    for(int t=0;t!=ntheta;t++)
    {
       for(int f=0;f!=nphi;f++)
       {
          sphere_dist[0][t*nphi+f]=acos(-1)*t/ntheta;
          sphere_dist[1][t*nphi+f]=2*acos(-1)*f/nphi;
       }
    }
       ofstream PICE;
//       ifstream sphere_dist_file;
//       sphere_dist_file.open("sphere_dist_file.txt");
/*       if(!sphere_dist_file.is_open())
       {
          std::cout<<"importing sphere distribution from previous comptutation"<<std::endl;
          ofstream write_sphere_dist_file;
          write_sphere_dist_file.open("sphere_dist_file.txt");
          for(int i=0;i!=n_points_sphere;i++)
          {
             std::cout<<"generating random point..."<<i<<"/"<<n_points_sphere<<std::endl;
             random=double(rand()%2000)-1000;
             random2=double(rand()%2000)-1000;
             random3=double(rand()%2000)-1000;
             temp=sqrt(pow(random,2)+pow(random2,2)+pow(random3,2));
             random/=temp;
             random2/=temp;
             random3/=temp;
             sphere_dist[0][i]=acos(random3);
             if(random>=0 && random2>=0)
               sphere_dist[1][i]=atan(random2/random);

             else if(random<0 && random2>=0)
               sphere_dist[1][i]=Pi+atan(random2/random);

             else if(random<0 && random2<0)
                sphere_dist[1][i]=Pi+atan(random2/random);

             else if(random>0 && random2 <0)
                sphere_dist[1][i]=2*Pi+atan(random2/random);
             write_sphere_dist_file<<sphere_dist[0][i]<<"  "<<sphere_dist[1][i]<<std::endl;
          }
          write_sphere_dist_file.close();
       }
       else
       {
          for(int i=0;i!=n_points_sphere;i++)
          {
            sphere_dist_file>>sphere_dist[0][i];
            sphere_dist_file>>sphere_dist[1][i];
          }
            sphere_dist_file.close();
       }*/
int m=11;
int n=0;
//for(int m=0;m!=n_states_neut;m++)
{
//   for(int n=0;n!=n_states_cat;n++)
   {
       ss_PICE.str("");
       ss_PICE<<"/data2/stephan/LiH_PICE_bruteforce_"<<m<<"_"<<n<<".txt";
       s_PICE=ss_PICE.str();

       cube_header(dyson_mo_basis_coeff,n_occ,n_states_neut,n_states_cat,neut_mo_cube_array,dyson_cube_loc,m,n,num_of_nucl,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,mo_cube_array); 

       vdMul(nx*ny*nz,coord[0],mo_cube_array,moment_cube_array[0]);
       vdMul(nx*ny*nz,coord[1],mo_cube_array,moment_cube_array[1]);
       vdMul(nx*ny*nz,coord[2],mo_cube_array,moment_cube_array[2]);
       std::cout<<"DYSON ORBITAL CUBE CONSTRUCTED"<<std::endl;

       PICE.open(s_PICE.c_str());

       for(int e=1;e!=nk+1;e++)
       {
          kp=e*kmax/nk;
  //        kp=0.3;
          sum=0;
//   std::cout<<"The factor is "<<4*0.441*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)*kp/(2*Pi*137)<<std::endl;

          for(int i=0;i!=n_points_sphere/continuum_matrix_cols;i++)
          {
             std::cout<<"Entering loop"<<std::endl;
             for(int j=0;j!=continuum_matrix_cols;j++)
             {
                theta[j]=sphere_dist[0][i*continuum_matrix_cols+j];
                phi[j]=sphere_dist[1][i*continuum_matrix_cols+j];
             }
             std::cout<<"building pw..."<<std::endl;

             pw_builder(Reinput,Iminput,kp,theta,phi,continuum_matrix_cols,xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz,n_occ,neut_mo_cube_array);
             std::cout<<"dipole moment along X...";
             cube_dot_product(Reinput,moment_cube_array[0],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Reproduct_vector[0]);
             cube_dot_product(Iminput,moment_cube_array[0],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Improduct_vector[0]);
             std::cout<<"dipole moment along Y...";
             cube_dot_product(Reinput,moment_cube_array[1],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Reproduct_vector[1]);
             cube_dot_product(Iminput,moment_cube_array[1],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Improduct_vector[1]);
             std::cout<<"dipole moment along Z...";
             cube_dot_product(Reinput,moment_cube_array[2],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Reproduct_vector[2]);
             cube_dot_product(Iminput,moment_cube_array[2],nx,ny,nz,(xmax-xmin)/nx,(ymax-ymin)/ny,(zmax-zmin)/nz,continuum_matrix_cols,Improduct_vector[2]);
          for(int w=0;w!=continuum_matrix_cols;w++)
          {
             PICE<<kp<<"    "<<theta[w]<<"    "<<phi[w]<<"    "<<Reproduct_vector[0][w]<<"   "<<Improduct_vector[0][w]<<"  "<<Reproduct_vector[1][w]<<"  "<<Improduct_vector[1][w]<<"   "<<Reproduct_vector[2][w]<<"  "<<Improduct_vector[2][w]<<std::endl;
             if((i*continuum_matrix_cols+w+1)%nphi==0)
             {
               PICE<<std::endl;
             }
          }
          }
       }
       PICE.close();
   }
}

    return 0;
    
}
