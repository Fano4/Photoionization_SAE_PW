//testing routines for algebra.cpp
//
//

//////////////////////////////////////////
//
//////////////////////////////////////////
bool test_determinant()
{
   bool test1(0);
   bool test2(0);
   bool test3(0);
   bool test4(0);

   std::cout<<" Testing determinant...";
   //Test_case_1
   {
     //Environment
     int rank1(10);
     double *matrix1;
     matrix1=new double [rank1*rank1];
      std::cout<<"1";
      //test_case
      for(int i=0;i!=rank1;i++)
         matrix1[i*rank1+i]=1;
      //test
      if(determinant(matrix1,rank1)==1)
         test1=1;
      //Finalize
      delete [] matrix1;
   }
   //Test_case_2
   {
     //Environment
     int rank2(3);
     double *matrix2;
     matrix2=new double [rank2*rank2];
     std::cout<<"2";
     double temp[9]={4,2,6,2,3,1,3,1,3};
     for(int i=0;i!=rank2*rank2;i++)
        matrix2[i]=temp[i];
        
      //test
      if(determinant(matrix2,rank2)==-16)
         test2=1;
      //Finalize
      delete [] matrix2;
   }
   //Test_case_3
   {
     //Environment
     double* matrix3;
     int rank3(5);
      std::cout<<"3";
      //test_case
      matrix3=new double [rank3*rank3];
      double temp3[25]={8,3,3,3,5,4,10,10,4,9,7,9,9,10,4,9,3,9,4,3,8,8,9,1,4};
     for(int i=0;i!=rank3*rank3;i++)
        matrix3[i]=temp3[i];
      //test
      if(determinant(matrix3,rank3)==17646)
         test3=1;
      //Finalize
      delete [] matrix3;
   }

   //Test_case_4
   {
     //Environment
     double* matrix4;
     int rank4(5);
      std::cout<<"4";
      //test_case
      matrix4=new double [rank4*rank4];
      double temp4[25]={67.2555,3.13869,6.82828,0.8213,5.70102,38.9396,4.17780,1.96137,30.5563,5.05692,86.5873,2.39583,0.65090,2.48812,9.72952,55.5325,5.85502,7.09267,2.91001,4.38487,.313663,.537773,.853602,8.3335,0.0001};
     for(int i=0;i!=rank4*rank4;i++)
        matrix4[i]=temp4[i];
      //test
      if( fabs(determinant(matrix4,rank4)+5261.323790016979)<=1e-10)
         test4=1;
//      std::cout<<std::setprecision(15)<<" "<<fabs(determinant(matrix4,rank4)+5261.323790016979)<<std::endl;
      //Finalize
      delete [] matrix4;
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
         std::cout<<"Error 2..";
      if(!test3)
         std::cout<<"Error 3..";
      if(!test4)
         std::cout<<"Error 4..";

      std::cout<<std::endl;
      return 0;
   }

}
//////////////////////////////////////////
//
//////////////////////////////////////////
bool test_prime_decomposer()
{

   std::cout<<" Testing test_prime_decomposer...";

   bool test1(1);
   int n(0);
   int val(1);
   srand (time(0));
   //MAX_FACTORIAL_PRIME = 25 maximum number of prime numbers supported => Max prime value = 97 
   //Test_case_1 : Recovering the right number
   {
      //environment
      for(int i=0;i!=10;i++)
      {
         std::cout<<i;
         //Environment
         int* n_prime=new int[MAX_FACTORIAL_PRIME];
         //test_case
          n=(rand() % 100 + 1);
          prime_decomposer(n,n_prime);
          val=1;
          for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
             val*=pow(PRIME[i],n_prime[i]);

          //test
          if(n == val)
             test1*=1;
          else
          {
             test1*=0;
             std::cout<<std::endl<<"failed for value "<<n<<" val = "<<val<<std::endl;
             for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
                std::cout<<n_prime[i]<<",";
          }
          //Finalize
          delete [] n_prime;
      }
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
/*bool test_factorized_sum()
{

   std::cout<<" Testing factorized_sum...";

   bool test1(1);
   int n(0);
   int val(1);
   srand (time(0));
   //MAX_FACTORIAL_PRIME = 25 maximum number of prime numbers supported => Max prime value = 97 
   //Test_case_1 : Recovering the right number
   {
      //environment
      int* a1_prime=new int[MAX_FACTORIAL_PRIME];
      int* a2_prime=new int[MAX_FACTORIAL_PRIME];
      int* sum_prime=new int[MAX_FACTORIAL_PRIME];
      int a1(0);
      int a2(0);
      int sum(0);
      for(int i=0;i!=10;i++)
      {
         std::cout<<i;
         a1=( rand() % 50 + 1);
         a2=( rand() % 50 + 1);
         //test case
         fact_prime_decomposer(a1,a1_prime);
         fact_prime_decomposer(a2,a2_prime);
         factorized_sum(a1_prime,a2_prime,sum_prime);
         sum=1;
         for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
             sum*=pow(PRIME[i],sum_prime[i]);
      //test
         if(sum==a1+a2)
            test1*=1;
         else
          {
             test1*=0;
             std::cout<<std::endl<<"failed for value "<<a1<<"+"<<a2<<" ; "<<" sum = "<<sum<<std::endl;
//             for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
//                std::cout<<n_prime[i]<<",";
          }

      }
      //Finalize
      delete [] a1_prime;
      delete [] a2_prime;
      delete [] sum_prime;

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
   return 0;
}*/
//////////////////////////////////////////
//
//////////////////////////////////////////
bool test_fact_prime_decomposer()
{
   std::cout<<" Testing test_fact_prime_decomposer...";

   bool test1(1);
   //MAX_FACTORIAL_PRIME = 25 maximum number of prime numbers supported for the decomposition
   //MAX_N_FACTORIAL = 100 Maximum number whose factorial can be computed
   //
   //Test_case_1 : The exponential series
   {
      std::cout<<1;
      //environment
      const int n_max(100);
      const double val_e(2.718281828459045); //15 decimals
      double sum=0;
      double approx_e(0);
      double val(1);
      const double accuracy(1e-10);
         //Environment
      int* n_fac_prime=new int[MAX_FACTORIAL_PRIME];

      for(int n=1;n!=n_max+1;n++)
      {
         val=1;

         //test_case

          fact_prime_decomposer(n,n_fac_prime);

          for(int i=0;i!=MAX_FACTORIAL_PRIME;i++)
             val*=pow(PRIME[i],-n_fac_prime[i]);

          approx_e+=val;

      }
          //test
          if(fabs(val_e-approx_e)/val_e <= accuracy)
          {
             test1*=1;
             std::cout<<"...passed"<<std::endl;
          }
          else
          {
             test1*=0;
             std::cout<<std::endl<<"failed. val = "<<approx_e<<std::endl;
          }
      delete [] n_fac_prime;
      
      //We wanto to compute a sum of inverse of factorial that are decomposed in prime numbers 
      //We first need to factorize all the factorials and then to compute the difference of the remain
      //Eventually, we factorize the remain and multiply the prefactor by the remain.

   }
}
//////////////////////////////////////////
//
//////////////////////////////////////////
/*
bool test_wigner3j()
{
   int l1,l2,l3,m1,m2,m3;
   bool test1(0);
   std::cout<<" Testing wigner3j...";

   srand (time(0));

   //Test_case_1
   {
      for(int i=0;i!=10;i++)
      {
         std::cout<<i;
         //Environment
          l1=(rand() % 10);
          l2=(rand() % 10);
          if( (l2+l1-fabs(l2-l1)) == 0)
             l3=l2+l1;
          else
             l3 = rand()% ( (l2+l1 ) - int(fabs(l2-l1))) + int(fabs(l2-l1)) ;
          m1=pow(-1,rand() % 2)*(rand() % (l1 + 1));
          m2=pow(-1,rand() % 2)*(rand() % (l2 + 1));
          m3=pow(-1,rand() % 2)*(rand() % (l3 + 1));
          std::cout<<std::endl<<" TEST: l1 = "<<l1<<", l2 = "<<l2<<", l3 = "<<l3<<", m1 = "<<m1<<", m2 = "<<m2<<", m3 = "<<m3<<std::endl;
          std::cout<<WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3)<<std::endl;

      }
      exit(EXIT_SUCCESS);
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
*/
