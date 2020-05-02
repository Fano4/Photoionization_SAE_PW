//testing routines for Legendre_functions.cpp
//
//
//////////////////////////////////////////
//
//////////////////////////////////////////
bool test_legendre()
{
   bool test1(0);
   bool test2(0);
   bool test3(0);

   double x[100]={-1.,-0.97979798,-0.95959596,-0.93939394,-0.91919192,
          -0.8989899 ,-0.87878788,-0.85858586,-0.83838384,-0.81818182,
          -0.7979798 ,-0.77777778,-0.75757576,-0.73737374,-0.71717172,
          -0.6969697 ,-0.67676768,-0.65656566,-0.63636364,-0.61616162,
          -0.5959596 ,-0.57575758,-0.55555556,-0.53535354,-0.51515152,
          -0.49494949,-0.47474747,-0.45454545,-0.43434343,-0.41414141,
          -0.39393939,-0.37373737,-0.35353535,-0.33333333,-0.31313131,
          -0.29292929,-0.27272727,-0.25252525,-0.23232323,-0.21212121,
          -0.19191919,-0.17171717,-0.15151515,-0.13131313,-0.11111111,
          -0.09090909,-0.07070707,-0.05050505,-0.03030303,-0.01010101,
           0.01010101, 0.03030303, 0.05050505, 0.07070707, 0.09090909,
           0.11111111, 0.13131313, 0.15151515, 0.17171717, 0.19191919,
           0.21212121, 0.23232323, 0.25252525, 0.27272727, 0.29292929,
           0.31313131, 0.33333333, 0.35353535, 0.37373737, 0.39393939,
           0.41414141, 0.43434343, 0.45454545, 0.47474747, 0.49494949,
           0.51515152, 0.53535354, 0.55555556, 0.57575758, 0.5959596 ,
           0.61616162, 0.63636364, 0.65656566, 0.67676768, 0.6969697 ,
           0.71717172, 0.73737374, 0.75757576, 0.77777778, 0.7979798 ,
           0.81818182, 0.83838384, 0.85858586, 0.87878788, 0.8989899 ,
           0.91919192, 0.93939394, 0.95959596, 0.97979798, 1.        };

   std::cout<<" Testing legendre...";
   //Test_case_1
   {
     //Environment
     unsigned int l(2);
     double *X1=new double[100];
     double *Y1=new double[100];
     //test_case
     std::cout<<"1";
     for(int i=0;i!=100;i++)
     {
        X1[i]=0.5*(3*x[i]*x[i]-1);
        Y1[i]=legendre(l,x[i]);
     }
     double r1(correlation_coefficient(X1,Y1,100));

     //test
     if(fabs(r1)>=0.999999)
         test1=1;
      //Finalize
      delete [] X1;
      delete [] Y1;
   }
   //
   //Test_case_2
   {
     //Environment
     unsigned int l(5);
     double *X2=new double[100];
     double *Y2=new double[100];
     //test_case
     std::cout<<"2";
     for(int i=0;i!=100;i++)
     {
        X2[i]=(1./8.)*(63*pow(x[i],5)-70*pow(x[i],3)+15*x[i]);
        Y2[i]=legendre(l,x[i]);
     }
     double r2(correlation_coefficient(X2,Y2,100));

     //test
     if(fabs(r2)>=0.999999)
         test2=1;
      //Finalize
      delete [] X2;
      delete [] Y2;
   }
   //Test_case_3
   {
     //Environment
     unsigned int l(8);
     double *X3=new double[100];
     double *Y3=new double[100];
     //test_case
     std::cout<<"3";
     for(int i=0;i!=100;i++)
     {
        X3[i]=(1./128.)*(6435*pow(x[i],8)-12012*pow(x[i],6)+6930*pow(x[i],4)-1260*pow(x[i],2)+35);
        Y3[i]=legendre(l,x[i]);
     }
     double r3(correlation_coefficient(X3,Y3,100));

     //test
     if(fabs(r3)>=0.999999)
         test3=1;
      //Finalize
      delete [] X3;
      delete [] Y3;
   }
   //Print results 
   if(test1*test2*test3)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      if(!test2)
         std::cout<<"Error 2...";
      if(!test3)
         std::cout<<"Error 3...";
      std::cout<<std::endl;
      return 0;
   }

}
//////////////////////////////////////////
//
//
//////////////////////////////////////////
bool test_associated_legendre_nonorm()
{
   bool test1(0);
   bool test2(0);
   bool test3(0);
   bool test4(0);

   double x[100]={-1.,-0.97979798,-0.95959596,-0.93939394,-0.91919192,
          -0.8989899 ,-0.87878788,-0.85858586,-0.83838384,-0.81818182,
          -0.7979798 ,-0.77777778,-0.75757576,-0.73737374,-0.71717172,
          -0.6969697 ,-0.67676768,-0.65656566,-0.63636364,-0.61616162,
          -0.5959596 ,-0.57575758,-0.55555556,-0.53535354,-0.51515152,
          -0.49494949,-0.47474747,-0.45454545,-0.43434343,-0.41414141,
          -0.39393939,-0.37373737,-0.35353535,-0.33333333,-0.31313131,
          -0.29292929,-0.27272727,-0.25252525,-0.23232323,-0.21212121,
          -0.19191919,-0.17171717,-0.15151515,-0.13131313,-0.11111111,
          -0.09090909,-0.07070707,-0.05050505,-0.03030303,-0.01010101,
           0.01010101, 0.03030303, 0.05050505, 0.07070707, 0.09090909,
           0.11111111, 0.13131313, 0.15151515, 0.17171717, 0.19191919,
           0.21212121, 0.23232323, 0.25252525, 0.27272727, 0.29292929,
           0.31313131, 0.33333333, 0.35353535, 0.37373737, 0.39393939,
           0.41414141, 0.43434343, 0.45454545, 0.47474747, 0.49494949,
           0.51515152, 0.53535354, 0.55555556, 0.57575758, 0.5959596 ,
           0.61616162, 0.63636364, 0.65656566, 0.67676768, 0.6969697 ,
           0.71717172, 0.73737374, 0.75757576, 0.77777778, 0.7979798 ,
           0.81818182, 0.83838384, 0.85858586, 0.87878788, 0.8989899 ,
           0.91919192, 0.93939394, 0.95959596, 0.97979798, 1.        };
   int lm[100][2]={{0,0},{1,-1},{1,0},{1,1},{2,-2},{2,-1},{2,0},{2,1},{2,2},{3,-3},{3,-2},{3,-1},{3,0},{3,1},{3,2},{3,3},{4,-4},{4,-3},{4,-2},{4,-1},{4,0},{4,1},{4,2},{4,3},{4,4},{5,-5},{5,-4},{5,-3},{5,-2},{5,-1},{5,0},{5,1},{5,2},{5,3},{5,4},{5,5},{6,-6},{6,-5},{6,-4},{6,-3},{6,-2},{6,-1},{6,0},{6,1},{6,2},{6,3},{6,4},{6,5},{6,6},{7,-7},{7,-6},{7,-5},{7,-4},{7,-3},{7,-2},{7,-1},{7,0},{7,1},{7,2},{7,3},{7,4},{7,5},{7,6},{7,7},{8,-8},{8,-7},{8,-6},{8,-5},{8,-4},{8,-3},{8,-2},{8,-1},{8,0},{8,1},{8,2},{8,3},{8,4},{8,5},{8,6},{8,7},{8,8},{9,-9},{9,-8},{9,-7},{9,-6},{9,-5},{9,-4},{9,-3},{9,-2},{9,-1},{9,0},{9,1},{9,2},{9,3},{9,4},{9,5},{9,6},{9,7},{9,8},{9,9}};

   std::cout<<" Testing associated_legendre_nonorm...";
   //Test_case_1
   {
     // l=1 m=-1
     //Environment
     unsigned int l(lm[1][0]);
     int m(lm[1][1]);
     double *X1=new double[100];
     double *Y1=new double[100];
     double *D1=new double[100];
     //test_case
     std::cout<<"1";
     for(int i=0;i!=100;i++)
     {
        X1[i]=0.5*pow(1-x[i]*x[i],0.5);
        Y1[i]=associated_legendre_nonorm(l,m,x[i]);
        D1[i]=(Y1[i]-X1[i]);
     }
     //double r1(correlation_coefficient(X1,Y1,100));
     double r1(1-sigma(D1,100));

     //test
     if(fabs(r1)>=0.999999)
         test1=1;
      //Finalize
      delete [] X1;
      delete [] Y1;
      delete [] D1;
   }
   //
   //Test_case_2
   {
     // l=1 m=1
     //Environment
     unsigned int l(lm[3][0]);
     int m(lm[3][1]);
     double *X2=new double[100];
     double *Y2=new double[100];
     double *D2=new double[100];
     //test_case
     std::cout<<"2";
     for(int i=0;i!=100;i++)
     {
        X2[i]=(pow(1-x[i]*x[i],0.5));
        Y2[i]=associated_legendre_nonorm(l,m,x[i]);
        D2[i]=Y2[i]-X2[i];
//        std::cout<<","<<D2[i]<<std::endl;
     }
     //double r2(correlation_coefficient(X2,Y2,100));
     double r2(1-sigma(D2,100));

     //test
     if(fabs(r2)>=0.999999)
         test2=1;
      //Finalize
      delete [] X2;
      delete [] Y2;
      delete [] D2;
   }
   //Test_case_3
   {
     // l=6 m=5
     //Environment
     unsigned int l(lm[47][0]);
     unsigned int m(lm[47][1]);
     double *X3=new double[100];
     double *Y3=new double[100];
     double *D3=new double[100];
     //test_case
     std::cout<<"3";
     for(int i=0;i!=100;i++)
     {
        X3[i]=(pow(1-x[i]*x[i],2.5)*x[i]*10395);
        Y3[i]=associated_legendre_nonorm(l,m,x[i]);
        D3[i]=Y3[i]-X3[i];
     }
//     double r3(correlation_coefficient(X3,Y3,100));
     double r3(1-sigma(D3,100));

     //test
     if(fabs(r3)>=0.999999)
         test3=1;
      //Finalize
      delete [] X3;
      delete [] Y3;
      delete [] D3;
   }
   //Test_case_4
   {
     // P_l^-m(x)=(-1)**m*(l-m)!/(l+m)! * P_l^m(x)
     //Environment
     unsigned int l(5);
     int m(5);
     double *X4=new double[100];
     double *Y4=new double[100];
     double *D4=new double[100];
     //test_case
     std::cout<<"4";
     for(int i=0;i!=100;i++)
     {
        X4[i]=associated_legendre_nonorm(l,m,x[i])*factorial(l-m)/factorial(l+m);
        Y4[i]=associated_legendre_nonorm(l,-m,x[i]);
        D4[i]=Y4[i]-X4[i];
     }
 //    double r4(correlation_coefficient(X4,Y4,100));
     double r4(1-sigma(D4,100));

     //test
     if(fabs(r4)>=0.999999)
         test4=1;
      //Finalize
      delete [] X4;
      delete [] Y4;
      delete [] D4;
   }
   //Print results 
   if(test1*test2*test3*test4)
   {
      std::cout<<"...passed"<<std::endl;
      return 1;
   }
   else
   {
      std::cout<<"...FAILED...";
      if(!test1)
         std::cout<<"Error 1...";
      if(!test2)
         std::cout<<"Error 2...";
      if(!test3)
         std::cout<<"Error 3...";
      std::cout<<std::endl;
      return 0;
   }

}
