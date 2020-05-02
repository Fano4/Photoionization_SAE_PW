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
