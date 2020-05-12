void build_transition_density_matrix(int n_states_neut,int n_closed,int n_occ,int ci_size_neut,int n_elec_neut,double **ci_vec_neut,double **tran_den_mat_mo)
{

   bool test(0);
   bool test2(0);
   bool test3(0);
   int q(0);
   int p(0);
   double sum(0);
   double det_val;
   
    for (int i=0; i!=n_states_neut; i++)//ELECTRONIC STATE N
    {
        for (int j=0; j!=n_states_neut; j++)//ELECTRONIC STATE K
        {
           std::cout<<" density between states "<<i<<" and "<<j<<std::endl;
           sum=0;
         for(int k=0;k!=(n_closed+n_occ);k++)
         {
            for(int kp=0;kp!=n_closed+n_occ;kp++)
            {
               tran_den_mat_mo[n_states_neut*i+j][(n_occ+n_closed)*k+kp] = 0; 
               for(int m=0;m!=ci_size_neut;m++)
               {
                  for(int n=0;n!=ci_size_neut;n++)
                  {

                     det_val=build_reduced_determinant(k,kp,n_elec_neut,n_closed,n_occ,&ci_vec_neut[0][(n_elec_neut+n_states_neut)*m],&ci_vec_neut[0][(n_elec_neut+n_states_neut)*n],&ci_vec_neut[1][n_elec_neut*m],&ci_vec_neut[1][n_elec_neut*n]);

                     tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]+=ci_vec_neut[0][(n_elec_neut+n_states_neut)*(m)+n_elec_neut+i]*ci_vec_neut[0][(n_elec_neut+n_states_neut)*(n)+n_elec_neut+j]*det_val;
                    }
                 }
               if(k==kp)
               {
                 sum+=tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp];
               //  std::cout<<" from orbital "<<k<<" and from orbital "<<kp<<":"<<tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]<<std::endl;
               }
//               std::cout<<std::setprecision(8)<<"trdm val "<<tran_den_mat_mo[i*n_states_neut+j][k*(n_occ+n_closed)+kp]<<std::endl;
              }
           }
         std::cout<<"SUM = "<<sum<<std::endl;
        }
    }
}
double build_reduced_determinant( int ai,int aj,int n_elec,int n_closed,int n_occ,double* mo_vector_1,double* mo_vector_2,double *spin_vector_1,double *spin_vector_2)
{
   /* Given the vectors containing the mo labels and the spin labels of the electrons, this routine builds a slater determinant from which one electron contained in the mo's i and j have been removed   !!!! ONLY FOR SINGLET AND SIMPLE EXCITATION
   */

   bool test2(0);
   bool test3(0);
   int spin(0);
   double temp(0);

   int new_vector_1[(n_occ+n_closed)];
   int new_vector_2[(n_occ+n_closed)];

   double prefactor(1);

   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      new_vector_1[k]=0;
      new_vector_2[k]=0;
   }
   
   for(int e=0;e!=n_elec;e++)
   {
      new_vector_1[int(mo_vector_1[e])]+=1;
      new_vector_2[int(mo_vector_2[e])]+=1;
   }
   /*
   std::cout<<"Taking electron from orbitals "<<ai<<","<<aj<<std::endl;
   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      std::cout<<new_vector_1[k]<<" ";
   }std::cout<<std::endl;
   
   for(int k=0;k!=(n_occ+n_closed);k++)
   {
      std::cout<<new_vector_2[k]<<" ";
   }std::cout<<std::endl;
   */
   prefactor=sqrt(double(new_vector_1[ai]))*sqrt(double(new_vector_2[aj]));
   new_vector_1[ai]--;
   new_vector_2[aj]--;

   for(int k=0;k!=(n_occ+n_closed);k++)
   {
//      std::cout<<new_vector_1[k]<<std::endl;
      if(new_vector_1[k]!=new_vector_2[k])
         return 0;
   }
   
  // if(prefactor != 0 )
  //    std::cout<<prefactor<<std::endl;
   return prefactor;
}
