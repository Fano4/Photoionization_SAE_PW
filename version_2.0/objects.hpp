#include <cmath>
class elec_struct {
   public:
      //Constructor
       elec_struct();
       elec_struct(int n_sym,int *n_states_neut,int *n_states_cat,int *n_occ,int *n_closed,int n_elec);
      //Set the parameters of the computation
      void set_n_sym(int n_sym);
      void set_n_states(int *n_states_neut,int *n_states_cat);
      void set_cas(int *n_occ,int *n_closed);
      void set_n_elec(int n_elec);
      void set_basis(int basis_size);
      void set_ci_size(int *ci_size_neut,int *ci_size_cat);
      void set_ci(int **ci_vec_neut,int **ci_vec_cat);
      void set_overlap(double *overlap);
      void set_mo_cube_loc(std::string cube_loc);
      //Return the parameters of the computation
      int  n_sym();
      int  n_states_neut(int sym);
      int  n_states_cat(int sym);
      int  n_closed(int sym);
      int  n_occ(int sym);
      int  n_elec();
      int  basis_size();
      int  ci_size_neut(int sym);
      int  ci_size_cat(int sym);
      int  ci_vector_neut(bool section, int index );
      int  ci_vector_cat(bool section, int index );
      double  overlap(int i,int j);

   private:
      //Member variables
      std::string m_mo_cube_loc;
      int m_n_sym;
      int *m_n_states_neut;
      int *m_n_states_cat;
      int *m_n_occ;
      int *m_n_closed;
      int m_n_elec_neut;
      int m_basis_size;
      int *m_ci_size_neut;
      int *m_ci_size_cat;
      double *m_ci_vec_neut[2];
      double *m_ci_vec_cat[2];
      double *m_overlap;
      //Private member functions

};
class cube_struct {
      //structure of the cubes and cube files used 
   public:
      //Constructor
      cube_struct();
      //set cube dimensions and properties
      void set_cube_dim(int nx,double xmin,double xmax,int ny, double ymin,double ymax,int nz,double zmin,double zmax);
      void set_titles(std::string pretitle,std::string title);
      void set_geom(int n_nucl,int *nucl_Z,double ** nucl_pos);
      void set_cube(double * cube);
      void make_header();

      //show cube dimensions and properties
      int nx();
      int ny();
      int nz();
      double xmin();
      double xmax();
      double ymin();
      double ymax();
      double zmin();
      double zmax();
      double cube(int i,int j,int k);

   private:

      std::string m_header;
      std::string *m_titles;

      int m_num_nucl;
      int *m_nucl_Z;
      double **m_nucl_pos;

      int m_nx;
      int m_ny;
      int m_nz;

      double m_xmin;
      double m_ymin;
      double m_zmin;
      double m_xmax;
      double m_ymax;
      double m_zmax;

      double *m_cube;

};
class continuum_struct {
   //structure of the continuum subspace
   public:
      //Constructor
      continuum_struct();

      void set_n_points_sphere(int n=0);
      
      void set_k_vector(int nk, int kmin, int kmax);
      void set_theta_vector(int ntheta);
      void set_phi_vector(int nphi);

      int nk();
      int ntheta();
      int nphi();
      int n_points_sphere();
      double kmin();
      double kmax();

   private:
      //dimension of the continuum (radius of the sphere)
      double m_kmin;
      double m_kmax;
      //number of points
      int m_nk;
      int m_ntheta;//theta=polar angle [0,PI]
      int m_nphi;//phi = azimuthal angle [0,2PI]
      int m_n_points_sphere; //number of points on the unit sphere = ntheta*nphi if the grid is structured ; can be chosen if the grid is randomly generated.

      //coordinates vectors
      double *m_modulus_k;
      double *m_theta;
      double *m_phi;
};
