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
    double **nucl_spher_pos;
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
    int ci_size_neut(0);
    int ci_size_cat(0);
    int ci_size_neut_sym[n_sym];
    int ci_size_cat_sym[n_sym];
    int n_elec_neut(4);//!!!! ecrire une routine qui cherche le nombre d'electrons dans l'output molpro!!!
   double Pi=acos(-1);
    double ***contraction_coeff_sym;
    double ***contraction_zeta_sym;
    int **contraction_number_sym;
    int **nucl_basis_func_sym;
    std::string** basis_func_type_sym;
    double **contraction_coeff;
    double **contraction_zeta;
    int *contraction_number;
    int *nucl_basis_func;
    std::string* basis_func_type;
    double* MO_coeff_neutral;
    const std::string molpro_output_path("molpro_output.out");


    //TEMPORARY VARIABLES
    double *overlap;
    double **mo_dipole;
    double **mo_dipole_spher;
    double *dyson_mo_basis_coeff;
    double *ci_vec_neut[2];
    double *ci_vec_cat[2];
    double sum(0);
    int index(0);
    int index2(0);
    int temp_int(0);



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
    std::cout<<"PROBE SIZE QUERY DONE"<<std::endl;//DEBOGAGE
   }
   else
      size_query(&n_occ,&n_closed,&basis_size,molpro_output_path);

    std::cout<<" number of occupied MO's : "<<n_occ<<std::endl;
    nucl_spher_pos=new double*[num_of_nucl];
    for(int i=0;i!=num_of_nucl;i++)
    {
       nucl_spher_pos[i]=new double[3];
    }
    nucl_spher_pos[0][0]=2.609900760;
    nucl_spher_pos[0][1]=acos(-1);
    nucl_spher_pos[0][2]=0.0;
    nucl_spher_pos[1][0]=0.437284517;
    nucl_spher_pos[1][1]=0.0;
    nucl_spher_pos[1][2]=0.0;


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
       std::cout<<"PROBE CI SIZE QUERY DONE"<<std::endl;//DEBOGAGE
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


    MO_coeff_neutral=new double[(n_occ+n_closed)*basis_size];
    overlap=new double[n_occ*n_occ];
    mo_dipole=new double *[3];
    mo_dipole_spher=new double *[3];
    mo_dipole[0]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole[1]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole[2]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[0]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[1]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
    mo_dipole_spher[2]=new double[(n_occ+n_closed)*(n_occ+n_closed)];
//*****************************COMPUTE DYSON ORBITALS*****************************
    //COMPUTE THE OVERLAP MATRIX BETWEEN THE MO OF THE NEUTRAL AND THE MO OF THE CATION FROM THE AO OVERLAP MATRIX
    if(symmetry)
       overlap_MO(overlap,n_occs,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral,n_sym);
    else
       overlap_MO(overlap,&n_occ,&basis_size,basis_size_sym,molpro_output_path,MO_coeff_neutral,n_sym);
    //COMPUTE THE MOLECULAR ORBITALS TRANSITION DIPOLE MOMENT MATRIX
    if(symmetry)
       dipole_MO(mo_dipole,n_occs,&basis_size,basis_size_sym,molpro_output_path,n_sym);
    else
       dipole_MO(mo_dipole,&n_occ,&basis_size,basis_size_sym,molpro_output_path,n_sym);

    for(int i=0;i!=n_occ*n_occ;i++)
    {
       mo_dipole_spher[0][i]=sqrt(mo_dipole[0][i]*mo_dipole[0][i]+mo_dipole[1][i]*mo_dipole[1][i]+mo_dipole[2][i]*mo_dipole[2][i]);
       if(mo_dipole_spher[0][i]==0)
       {
          mo_dipole_spher[1][i]=0;
          mo_dipole_spher[2][i]=0;
       }
       else if(mo_dipole[0][i]==0)
       {
          mo_dipole_spher[1][i]=acos(mo_dipole[2][i]/mo_dipole_spher[0][i]);
          if(mo_dipole[1][i]==0)
          {
             mo_dipole_spher[2][i]=0;
          }
          else if (mo_dipole[1][i]>0)
          {
             mo_dipole_spher[2][i]=acos(-1)/2;
          }
          else
          {
             mo_dipole_spher[2][i]=3*acos(-1)/2;
          }
       }
       else
       {
           mo_dipole_spher[1][i]=acos(mo_dipole[2][i]/mo_dipole_spher[0][i]);
           mo_dipole_spher[2][i]=atan(mo_dipole[1][i]/mo_dipole[0][i]);
       }
//       std::cout<<"["<<int(i/n_occ)<<","<<int(i%n_occ)<<"]"<<"("<<mo_dipole[0][i]<<","<<mo_dipole[1][i]<<","<<mo_dipole[2][i]<<") => "<<"("<<mo_dipole_spher[0][i]<<","<<mo_dipole_spher[1][i]<<","<<mo_dipole_spher[2][i]<<")"<<std::endl;
    }
    //GET THE INFORMATION ABOUT THE BASIS SET
    contraction_coeff_sym=new double**[n_sym];
    contraction_zeta_sym=new double**[n_sym];
    contraction_number_sym=new int*[n_sym];
    nucl_basis_func_sym=new int*[n_sym];
    basis_func_type_sym=new std::string*[n_sym];

    int total(0);

    for(int s=0;s!=n_sym;s++)
    {
       contraction_number_sym[s]=new int[basis_size_sym[s]];
       contraction_coeff_sym[s]=new double*[basis_size_sym[s]];
       contraction_zeta_sym[s]=new double*[basis_size_sym[s]];
       nucl_basis_func_sym[s]=new int[basis_size_sym[s]];
       basis_func_type_sym[s]=new std::string[basis_size_sym[s]];
       total+=basis_size_sym[s];

    }
    contraction_coeff=new double*[total];
    contraction_zeta=new double*[total];
    contraction_number=new int[total];
    nucl_basis_func=new int[total];
    basis_func_type=new std::string[total];

    basis_size_data_reader(n_sym, basis_size_sym,contraction_number_sym,molpro_output_path);
    total=0;
    for(int s=0;s!=n_sym;s++)
    {
       for(int t=0;t!=basis_size_sym[s];t++)
       {
          contraction_coeff[total]=new double[contraction_number_sym[s][t]];
          contraction_zeta[total]=new double[contraction_number_sym[s][t]];
          contraction_coeff_sym[s][t]=new double [contraction_number_sym[s][t]];
          contraction_zeta_sym[s][t]=new double [contraction_number_sym[s][t]];
//          std::cout<<contraction_number_sym[s][t]<<" , "<<s<<" , "<<t<<std::endl;
          total++;
       }
    }
    basis_data_reader(n_sym,basis_size_sym,contraction_number_sym,contraction_coeff_sym,contraction_zeta_sym,nucl_basis_func_sym,basis_func_type_sym,molpro_output_path);
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
//          std::cout<<basis_func_type[total]<<" !!"<<std::endl;
          total++;
       }
    }
    double **contraction_coeff_array=new double *[basis_size];
    double **contraction_zeta_array=new double *[basis_size];
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
//    exit(EXIT_SUCCESS);
/*
    for(int i=0;i!=n_occ;i++)
    {
       for(int j=0;j!=n_occ;j++)
       {
          std::cout<<overlap[i*n_occ+j]<<"   ";
       }std::cout<<std::endl;
    }*/
    //COMPUTE THE DYSON MO COEFFICIENTS IN THE BASIS OF THE MO OF THE NEUTRAL
    dyson_mo_coeff_comp( n_states_neut,n_states_cat, n_occ,ci_size_neut, ci_size_cat, n_elec_neut, ci_vec_neut, ci_vec_cat,overlap, dyson_mo_basis_coeff);

    for(int i=0;i!=n_states_neut;i++)
    {
       for(int j=0;j!=n_states_cat;j++)
       {
          temp_norm=0;
          for(int t=0;t!=n_occ;t++)
          {
             temp_norm+=dyson_mo_basis_coeff[i*n_states_cat*n_occ+j*n_occ+t]*dyson_mo_basis_coeff[i*n_states_cat*n_occ+j*n_occ+t];
          }
          std::cout<<"Dyson orbital norm between states "<<i<<" and "<<j<<" : "<<temp_norm<<std::endl;
       }
    }

    /*
     TESTING HF5 DIALOG
     */
      double nucl_coord(1.6125/0.529);
      int grid_size(1);
      int nucl_dim(1);
      write_output("test.h5", &n_states_neut, &n_states_cat, &n_occ, &n_closed,1,&grid_size,&nucl_coord,2,nucl_spher_pos,mo_dipole,&basis_size,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_number,contraction_coeff,contraction_zeta,nucl_basis_func,basis_func_type);
      read_output("test.h5", &n_states_neut, &n_states_cat, &n_occ, &n_closed,&nucl_dim,&grid_size,&nucl_coord,&num_of_nucl,nucl_spher_pos,mo_dipole,&basis_size,MO_coeff_neutral,dyson_mo_basis_coeff,contraction_number,contraction_coeff,contraction_zeta,nucl_basis_func,basis_func_type);
      exit(EXIT_SUCCESS);
    /*
     TESTING HF5 DIALOG
     */

    double kp(0);
    int state_neut(1);
    double int_cs(0);
    double thet(0);
    double phi(0);
    int nk=200;
    std::complex<double> temp;
    std::complex<double> sinet;
    std::complex<double> cosinet;
    std::complex<double> sinep;
    std::complex<double> cosinep;
    std::complex<double> modulus;
    int n_theta=50;
    int n_phi=40;
    //double efield[3]={1,0,0};//SPHERICAL COORDINATES
    double efield[3]={0,0,0.05};//CARTESIAN COORDINATES
    int i(0);
    int t(0);
    int p(0);
    int j(0);
//      MO_Fourier_transform(0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size);
/*
    for(int k=0;k!=nk;k++)
    {
       kp=2.7*(k+1)/nk;
       thet=0;
       phi=0;
       temp=0;
      for(int i=0;i!=n_occ;i++)
      {
          temp-=std::complex<double>(0,1)*dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
             *std::complex<double>(
             real(MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type, MO_coeff_neutral,basis_size))*efield[0]
                 *(sin(real(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)))*sin(efield[1])
                 *cos(real(MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))-efield[2])
                 +cos(real(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)))*cos(efield[1]))
                 ,imag(MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type, MO_coeff_neutral,basis_size))*efield[0]
                 *(sin(imag(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)))*sin(efield[1])
                    *cos(imag(MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))-efield[2])
                    +cos(imag(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)))*cos(efield[1])));
      }
       std::cout<<kp*kp*27.211/2<<","<<abs(temp)*0.05*0.05*27.211<<std::endl;
    }
*/

    for(int k=0;k!=nk;k++)
    {
       kp=1.35541*(k+1)/nk;
       int_cs=0;
#pragma omp parallel for reduction(+:int_cs) private (i,j,t,p,temp,thet,phi)
       for( t=0;t!=n_theta;t++)
       {
          thet=t*(acos(-1))/n_theta;
          for( p=0;p!=n_phi;p++)
          {
             phi=2*p*acos(-1)/n_phi;
             temp=0;
             for( i=0;i!=n_occ;i++)
             {
                 temp-=std::complex<double>(0,1)*dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
                          *(
                                efield[0]*(
                                   sin(thet)*cos(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   +cos(thet)*cos(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                -sin(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                                +efield[1]*(
                                   sin(thet)*sin(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   +cos(thet)*sin(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                -cos(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                                +efield[2]*(
                                   cos(thet)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   -sin(thet)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                           );
          
      /*          for(int j=0;j!=n_occ;j++)
                {
                
                   temp-=dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *MO_Fourier_transform(j,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)*(mo_dipole[0][i*n_occ+j]*efield[0]+mo_dipole[1][i*n_occ+j]*efield[1]+mo_dipole[2][i*n_occ+j]*efield[2]);
                 }*/
 //            std::cout<<temp<<std::endl;
             }
             

//             temp=MO_Fourier_transform(1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size);
             int_cs+=kp*sin(thet)*(acos(-1)/n_theta)*(2*acos(-1)/n_phi)*abs(temp)*abs(temp);
         //std::cout<<kp<<" , "<<thet<<" , "<<phi<<" , "<<temp<<std::endl;
//             std::cout<<"small loop "<<p<<std::endl;
          }
//             std::cout<<"LARGE loop "<<t<<std::endl<<"##################################################################"<<std::endl;
       }
       
          std::cout<<kp*kp*27.211/2<<","<<int_cs*27.211<<std::endl;
          
    }
    
/*

    std::cout<<efield[0]<<","<<efield[1]<<","<<efield[2]<<std::endl;
    kp=0.67;
       for(int t=0;t!=n_theta;t++)
       {
          thet=t*acos(-1)/n_theta;
          for(int p=0;p!=n_phi;p++)
          {
             phi=2*p*acos(-1)/n_phi;
             temp=0;
             //i=3;
//             std::cout<<thet<<"   "<<phi<<"   "<<real(MO_Fourier_transform(i,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))<<"   "<<imag(MO_Fourier_transform(i,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))<<std::endl;
             //std::cout<<thet<<"   "<<phi<<"   "<<real(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))<<"   "<<imag(MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))<<std::endl;
             //continue;
             for(int i=0;i!=n_occ;i++)
             {
                 temp-=std::complex<double>(0,1)*dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
                          *(
                                efield[0]*(
                                   sin(thet)*cos(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   +cos(thet)*cos(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                -sin(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                                +efield[1]*(
                                   sin(thet)*sin(phi)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   +cos(thet)*sin(phi)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                -cos(phi)*MO_Fourier_transform_grad(i,2,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                                +efield[2]*(
                                   cos(thet)*MO_Fourier_transform_grad(i,0,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)
                                   -sin(thet)*MO_Fourier_transform_grad(i,1,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size))
                           );
          
                for(int j=0;j!=n_occ;j++)
                {
                
                   temp-=dyson_mo_basis_coeff[state_neut*n_states_cat*n_occ+0*n_occ+i]
                      *MO_Fourier_transform(j,kp,thet,phi,nucl_spher_pos,nucl_basis_func,contraction_number,contraction_coeff,contraction_zeta,basis_func_type,MO_coeff_neutral,basis_size)*(mo_dipole[0][i*n_occ+j]*efield[0]+mo_dipole[1][i*n_occ+j]*efield[1]+mo_dipole[2][i*n_occ+j]*efield[2]);
                 }
             }
             std::cout<<thet<<"  "<<phi<<"   "<<real(temp)<<"    "<<imag(temp)<<std::endl;
          }std::cout<<std::endl;
       }
  */     
/*    std::cout<<"DYSON COEFF ROUTINE ENDED WITHOUT ISSUE"<<std::endl;

    double dtemp(0);
    for(int a=0;a!=n_states_neut;a++)
    {
       dtemp=0;
       for(int k=0;k!=n_occ;k++)
       {
          dtemp+=dyson_mo_basis_coeff[a*n_states_cat*n_occ+0*n_occ+k]*dyson_mo_basis_coeff[a*n_states_cat*n_occ+0*n_occ+k];
       }
       std::cout<<"Norm of Dyson orbital "<<a<<"_"<<0<<" = "<<dtemp<<std::endl;
       */
//cube_header(dyson_mo_basis_coeff,n_occ,n_states_neut,n_states_cat,neut_mo_cube_array,dyson_cube_loc.c_str(),a,0,2,nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,mo_cube_array);
      return 0;
    }

