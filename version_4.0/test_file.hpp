#ifndef test_file_hpp
#define test_file_hpp

double test2_integral(int l1,int l2,int l3,int m1,int m2,int m3,double* lnfact_memo);
void test_radial(int l1,int l2,int l3,int m1,int m2,int m3,double zeta,double kp,double* r0,double* lnfact_memo);
void cart_to_spher(double* x,double* y,double* z,double * r,double* t,double *f);
double numerical_integral(int l1,int m1,int l2,int m2,double zeta,double kp,double* r0,double* lnfact_memo);



#endif
