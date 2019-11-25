#ifndef test_file_hpp
#define test_file_hpp

double test2_integral(int l1,int l2,int l3,int m1,int m2,int m3);
void test_radial(int l1,int l2,int l3,int m1,int m2,int m3,double zeta,double kp,double* r0);
void cart_to_spher(double* x,double* y,double* z,double * r,double* t,double *f);
double numerical_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0);
void pw_bessel_overlap_comparison(int l2,int m2,double zeta,double kp,double thet,double phi,double* r0);
void pw_bessel_comparison(double kp,double kthet,double kphi,double r,double thet,double phi);



#endif
