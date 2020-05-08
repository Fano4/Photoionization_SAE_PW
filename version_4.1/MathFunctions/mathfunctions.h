double azim_integ(int m1);
double two_azim_integ(int m1,int m2);
double three_azim_integ(int m1,int m2,int m3);

void Jint_sort_indices(int* l1,int* l2,int* l3,int* m1,int* m2,int* m3);
double Jint_signflip_renormalize(int l1,int l2,int l3,int* m1,int* m2,int* m3);
double ALP_normalize(int l,int m);
double Jint_normalize(int l1,int l2,int l3,int m1,int m2,int m3);
bool Jint_special_cases(int l1,int l2,int l3,int m1,int m2,int m3,double* result);
double ALP_integral(int l,int m);
double two_ALP_integral(int l1,int l2,int m1,int m2);
double three_ALP_J_integral(int l1,int l2,int l3,int m1,int m2,int m3);

double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_p1(int l1,int l2,int l3,int m1,int m2,int m3);
double J_int_m1_D(int l1,int l2,int l3,int m1,int m2,int m3);

double I_m1_integral(int m1,int m2,int m3);
double I_p1_integral(int m1,int m2,int m3);
double I_m1_D_integral(int m1,int m2,int m3);
double I_p1_D_integral(int m1,int m2,int m3);

