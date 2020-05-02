//testing routines for angular_int_aux.cpp
//
//

//////////////////////////////////////////
//
//////////////////////////////////////////
/*bool azim_integ_test()
{
   int m1,m2,m3;

   srand (time(0));
   m1=pow(-1,rand() % 2)*(rand() % 10);

   for(int i=0;i!=100;i++)
   {
       m1=pow(-1,rand() % 2)*(rand() % 10);
       m2=pow(-1,rand() % 2)*(rand() % 10);
       m3=pow(-1,rand() % 2)*(rand() % 10);

       std::cout<<m1<<","<<m2<<","<<m3<<std::endl;
   }
   exit(EXIT_SUCCESS);
   return 0;
}*/
//////////////////////////////////////////
//
//double J_int_m2(int l1,int l2,int l3,int m1,int m2,int m3)
//
//////////////////////////////////////////
bool J_int_m2_test()
{
   bool test1(1);
   double threshold(1e-7);
   std::cout<<" Testing J_int_m2_test...";

   //Test_case_1
   {
      int nx(20000);
      double dx(2./nx);
      double x;
      srand (time(0));

      int l1,l2,l3,m1,m2,m3,nz1,nz2,nz3,nzeros;
      double sum(0);

      for(int i=0;i!=100;i++)
      {
         std::cout<<i;
         //Environment
          l1=(rand() % 10);
          l2=(rand() % 10);
          l3=(rand() % 10);
          m1=pow(-1,rand() % 2)*(rand() % (l1 + 1));
          m2=pow(-1,rand() % 2)*(rand() % (l2 + 1));
          m3=pow(-1,rand() % 2)*(rand() % (l3 + 1));
          std::cout<<std::endl<<" TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;

          nz1=(l1-fabs(m1));
          nz2=(l2-fabs(m2));
          nz3=(l3-fabs(m3));

          nzeros=pow(nz1,bool(nz1!=0))*pow(nz2,bool(nz2!=0))*pow(nz3,bool(nz3!=0));

          std::cout<<"NNNN"<<nzeros<<std::endl;
          //test_case
          sum=0;
          for(int n=0;n!=nx;n++)
          {
             x=-0.9999999999+n*1.9999999999/nx;
             sum+=dx*((1/sqrt(1-x*x))*associated_legendre_nonorm(l1,m1,x)*associated_legendre_nonorm(l2,m2,x)*associated_legendre_nonorm(l3,m3,x));
          }
          //test
          if(fabs(J_int_m2(l1,l2,l3,m1,m2,m3)-sum) <= 1e-10)
             test1*=1;
          else 
          {
             test1*=0;
             std::cout<<std::endl<<"FAILED TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;
             std::cout<<J_int_m2(l1,l2,l3,m1,m2,m3)<<" - "<<sum<<" >1e-10"<<std::endl;
          }
      }
      
      //Finalize
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
//////////////////////////////////////////
//
//////////////////////////////////////////
bool J_int_m1_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool J_int_p1_test(); 
//////////////////////////////////////////
//
//////////////////////////////////////////
bool J_int_D_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool J_int_m1_D_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool J_int_p1_D_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool I_m1_integral_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool I_p1_integral_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool I_m1_D_integral_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool I_p1_D_integral_test();
//////////////////////////////////////////
//
//////////////////////////////////////////
bool gaunt_formula_test()
{
   bool test1(1);
   std::cout<<" Testing gaunt_formula...";

   //Test_case_1
   {
      int nx(200000);
      double dx(2./nx);
      double threshold(dx);
      double x;
      srand (time(0));

      int l1,l2,l3,m1,m2,m3,nz1,nz2,nz3,nzeros;
      double sum(0);

      for(int i=0;i!=10;i++)
      {
         std::cout<<i;
         //Environment
          l1=(rand() % 10);
          l2=(rand() % 10);
          l3=(rand() % 10);
          m1=pow(-1,rand() % 2)*(rand() % (l1 + 1));
          m2=pow(-1,rand() % 2)*(rand() % (l2 + 1));
          m3=pow(-1,rand() % 2)*(rand() % (l3 + 1));
//          std::cout<<std::endl<<" TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;

//          nz1=(l1-fabs(m1));
//          nz2=(l2-fabs(m2));
//          nz3=(l3-fabs(m3));

//          nzeros=pow(nz1,bool(nz1!=0))*pow(nz2,bool(nz2!=0))*pow(nz3,bool(nz3!=0));

//          std::cout<<"NNNN"<<nzeros<<std::endl;
          //test_case
          sum=0;
          for(int n=0;n!=nx;n++)
          {
             x=-0.9999999999+n*1.9999999999/nx;
             sum+=dx*(associated_legendre_nonorm(l1,m1,x)*associated_legendre_nonorm(l2,m2,x)*associated_legendre_nonorm(l3,m3,x));
          }
          //test
          //
          if(gaunt_formula(l1,l2,l3,m1,m2,m3)!=0)
          {
             if(fabs(gaunt_formula(l1,l2,l3,m1,m2,m3)-sum)/fabs(gaunt_formula(l1,l2,l3,m1,m2,m3)) <= 1e-10)
                test1*=1;
             else 
             {
                test1*=0;
                std::cout<<std::endl<<"FAILED TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;
                std::cout<<gaunt_formula(l1,l2,l3,m1,m2,m3)<<" - "<<sum<<std::endl;
             }
          }
          else
          {
             if(fabs(sum)<=dx)
                test1*=1;
             else
             {
                test1*=0;
                std::cout<<std::endl<<"FAILED TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;
                std::cout<<gaunt_formula(l1,l2,l3,m1,m2,m3)<<" - "<<sum<<std::endl;
             }
          }
      }
      
      //Finalize
   }
   if(test1)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      std::cout<<std::endl;
      return 0;
   }
}
