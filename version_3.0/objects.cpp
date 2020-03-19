//elec struct class member functions
//
//
//
elec_struct::elec_struct()
{
   int temp(0);
   this->set_n_sym(1);
   this->set_n_states(&temp,&temp);
   this->set_cas(&temp,&temp);
   this->set_n_elec(0);

}
void elec_struct::set_n_sym(int n_sym)
{
   this->m_n_sym=n_sym;
}
void elec_struct::set_n_states(int *n_states_neut,int *n_states_cat)
{
   this->m_n_states_neut=new int [m_n_sym];
   this->m_n_states_cat=new int [m_n_sym];
   for(int i=0;i!=m_n_sym;i++)
   {
      this->m_n_states_neut[i]=n_states_neut[i];
      this->m_n_states_cat[i]=n_states_cat[i];
   }
}
void elec_struct::set_cas(int *n_occ,int *n_closed)
{
   this->m_n_occ=new int[m_n_sym];
   this->m_n_closed=new int[m_n_sym];
   for(int i=0;i!=m_n_sym;i++)
   {
      this->m_n_occ[i]=n_occ[i];
      this->m_n_closed[i]=n_closed[i];
   }

}
void elec_struct::set_n_elec(int n_elec)
{
   this->m_n_elec_neut=n_elec;
}
void elec_struct::set_basis(int basis_size)
{
   this->m_basis_size=basis_size;
}
void elec_struct::set_ci_size(int *ci_size_neut,int *ci_size_cat)
{
   this->m_ci_size_neut=new int[this->m_n_sym];
   this->m_ci_size_cat=new int[this->m_n_sym];
   for(int i=0;i!=this->m_n_sym;i++)
   {
      this->m_ci_size_neut[i]=ci_size_neut[i];
      this->m_ci_size_cat[i]=ci_size_cat[i];
   }
}
void elec_struct::set_ci(int **ci_vec_neut,int **ci_vec_cat)
{
   int tot_ci_size_neut(0);
   int tot_ci_size_cat(0);
   int tot_n_states_neut(0);
   int tot_n_states_cat(0);
   for(int i=0;i!=this->m_n_sym;i++)
   {
      tot_ci_size_neut+=this->m_ci_size_neut[i];
      tot_ci_size_cat+=this->m_ci_size_cat[i];
      tot_n_states_neut+=this->m_n_states_neut[i];
      tot_n_states_cat+=this->m_n_states_cat[i];
   }
   this->m_ci_vec_neut[0]=new double[this->m_n_elec_neut*tot_ci_size_neut+tot_n_states_neut*tot_ci_size_neut];
   this->m_ci_vec_neut[1]=new double[this->m_n_elec_neut*tot_ci_size_neut];
   this->m_ci_vec_cat[0]=new double[(this->m_n_elec_neut-1)*tot_ci_size_cat+tot_n_states_cat*tot_ci_size_cat];
   this->m_ci_vec_cat[1]=new double[(this->m_n_elec_neut-1)*tot_ci_size_cat];

   for(int i=0;i!=m_n_elec_neut*tot_ci_size_neut+tot_n_states_neut*tot_ci_size_neut;i++)
   {
      this->m_ci_vec_neut[0][i]=ci_vec_neut[0][i];
   }
   for(int i=0;i!=(m_n_elec_neut-1)*tot_ci_size_cat+tot_n_states_cat*tot_ci_size_cat;i++)
   {
      this->m_ci_vec_cat[0][i]=ci_vec_cat[0][i];
   }
   for(int i=0;i!=m_n_elec_neut*tot_ci_size_neut;i++)
   {
      this->m_ci_vec_neut[1][i]=ci_vec_neut[1][i];
   }
   for(int i=0;i!=(m_n_elec_neut-1)*tot_ci_size_cat;i++)
   {
      this->m_ci_vec_cat[1][i]=ci_vec_cat[1][i];
   }
}
void elec_struct::set_overlap(double *overlap)
{
   this->m_overlap=new double[this->m_basis_size*this->m_basis_size];
   for(int i=0;i!=this->m_basis_size*this->m_basis_size;i++)
   {
      this->m_overlap[i]=overlap[i];
   }
}
int elec_struct::n_sym()
{
   return this->m_n_sym;
}
int elec_struct::n_states_neut(int sym)
{
   return this->m_n_states_neut[sym];
}
int elec_struct::n_states_cat(int sym)
{
   return this->m_n_states_cat[sym];
}
int  elec_struct::n_closed(int sym)
{
   return this->m_n_closed[sym];
}
int  elec_struct::n_occ(int sym)
{
   return this->m_n_occ[sym];
}
int  elec_struct::n_elec()
{
   return this->m_n_elec_neut;
}
int  elec_struct::basis_size()
{
   return this->m_basis_size;
}
int  elec_struct::ci_size_neut(int sym)
{
   return this->m_ci_size_neut[sym]; 
}
int  elec_struct::ci_size_cat(int sym)
{
   return this->m_ci_size_cat[sym];
}
int  elec_struct::ci_vector_neut(bool section, int index )
{
   return this->m_ci_vec_neut[section][index];
}
int  elec_struct::ci_vector_cat(bool section, int index )
{
   return this->m_ci_vec_cat[section][index];
}
double elec_struct::overlap(int i,int j)
{
   return this->m_overlap[i*m_basis_size+j];
}
//
//
//continuum_struct class member functions
//
//
//
int continuum_struct::nk()
{
   return this->m_nk;
}
int continuum_struct::ntheta()
{
   return this->m_ntheta;
}
int continuum_struct::nphi()
{
   return this->m_nphi;
}
int continuum_struct::n_points_sphere()
{
   return this->m_n_points_sphere;
}
double continuum_struct::kmin()
{
   return  this->m_kmin;
}
double continuum_struct::kmax()
{
   return this->m_kmax;
}
void continuum_struct::set_k_vector(int nk, int kmin, int kmax)
{
   double k;
   this->m_modulus_k=new double [nk];
   for(int i=0;i!=nk;i++)
   {
      k=kmin+i*(kmax-kmin/nk);
      this->m_modulus_k[i]=k;
   }
   this->m_kmin=kmin;
   this->m_kmax=kmax;
   this->m_nk=nk;
}
void continuum_struct::set_theta_vector(int ntheta)
{
   double theta;
   this->m_theta=new double [ntheta];
   for(int i=0;i!=ntheta;i++)
   {
      theta=i*acos(-1)/ntheta;
      this->m_theta[i]=theta;
   }
   this->m_ntheta=ntheta;
}
void continuum_struct::set_phi_vector(int nphi)
{
   double phi;
   this->m_phi=new double [nphi];
   for(int i=0;i!=nphi;i++)
   {
      phi=i*2*acos(-1)/nphi;
      this->m_phi[i]=phi;
   }
   this->m_nphi=nphi;
}
void continuum_struct::set_n_points_sphere(int n)
{
   this->m_n_points_sphere=n;
}
//
//
//
//cube_struct class member functions
//
//
//
void cube_struct::set_cube_dim(int nx,double xmin,double xmax,int ny, double ymin,double ymax,int nz,double zmin,double zmax)
{
   this->m_nx=nx;
   this->m_ny=ny;
   this->m_nz=nz;
   this->m_xmin=xmin;
   this->m_xmax=xmax;
   this->m_ymin=ymin;
   this->m_ymax=ymax;
   this->m_zmin=zmin;
   this->m_zmax=zmax;
}
void cube_struct::set_titles(std::string pretitle,std::string title)
{
   this->m_titles[0]=pretitle;
   this->m_titles[1]=title;
}
void cube_struct::set_geom(int n_nucl,int *nucl_Z,double ** nucl_pos)
{
   this->m_num_nucl=n_nucl;
   this->m_nucl_Z=new int [n_nucl];
   this->m_nucl_pos=new double * [n_nucl];
   for(int i=0;i!=n_nucl;i++)
   {
      this->m_nucl_pos[i]=new double [3];
      for(int j=0;j!=2;j++)
      {
         this->m_nucl_pos[i][j]= nucl_pos[i][j];
      }
   }
}
void cube_struct::set_cube(double * cube)
{
   this->m_cube=new double[m_nx*m_ny*m_nz];
   for(int i=0;i!=m_nx;i++)
   {
      for(int j=0;j!=m_ny;j++)
      {
         for(int k=0;k!=m_nz;k++)
         {
            this->m_cube[i*m_ny*m_nz+j*m_nz+k]=cube[i*m_ny*m_nz+j*m_nz+k];
         }
      }
   }
}
int cube_struct::nx()
{
   return this->m_nx;
}
int cube_struct::ny()
{
   return this->m_ny;
}
int cube_struct::nz()
{
   return this->m_nz;
}
double cube_struct::xmin()
{
   return this->m_xmin;
}
double cube_struct::ymin()
{
   return this->m_ymin;
}
double cube_struct::zmin()
{
   return this->m_zmin;
}
double cube_struct::xmax()
{
   return this->m_xmax;
}
double cube_struct::ymax()
{
   return this->m_ymax;
}
double cube_struct::zmax()
{
   return this->m_zmax;
}
double cube_struct::cube(int i,int j,int k)
{
   return this->m_cube[i*m_ny*m_nz+j*m_nz+k];
}
void cube_struct::make_header()
{
   using namespace std;
   stringstream stream;
   stream.str("");
   stream<<this->m_titles[0]<<endl<<this->m_titles[1]<<endl;
   stream<<"-"<<this->m_num_nucl<<"  "<<this->m_xmin<<"  "<<this->m_ymin<<"  "<<this->m_zmin<<endl;
   stream<<this->m_nx<<"  "<<(this->m_xmax-this->m_xmin)/this->m_nx<<"  "<<"0.000000"<<"  "<<"0.000000"<<endl;
   stream<<this->m_ny<<"  "<<"0.000000"<<"  "<<(this->m_ymax-this->m_ymin)/this->m_ny<<"  "<<"0.000000"<<endl;
   stream<<this->m_nz<<"  "<<"0.000000"<<"  "<<"0.000000"<<"  "<<(this->m_zmax-this->m_zmin)/this->m_nz<<endl;
   for(int i=0;i!=this->m_num_nucl;i++)
   {
     stream<<this->m_nucl_Z[i]<<"   "<<double(this->m_nucl_Z[i]);
     for(int j=0;j!=3;j++)
     {
        stream<<this->m_nucl_pos[i][j]<<"   ";
     }stream<<endl;
   }
   stream<<"1   111"<<endl;
   this->m_header=stream.str();
}
