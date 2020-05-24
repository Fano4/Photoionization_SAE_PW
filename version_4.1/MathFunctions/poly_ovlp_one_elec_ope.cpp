#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "utilities.h"
#include "mathfunctions.h"

double slater_ovlp(std::vector<int> mo_a,std::vector<int> mo_b,std::vector<int> spin_a,std::vector<int> spin_b,std::vector<double> MO_S)
{
   //This function computes the overlap between two slater determinants without computing a determinant. 
   //This version requires the two determinants to have the same size
   
   using namespace std;

   if(mo_a.size()!=mo_b.size())
      err_diff_on_vec_size(mo_a.size(),mo_b.size());

   double sum(0);
   int n_occ(sqrt(MO_S.size()));
   
   double* matrix=new double[mo_a.size()*mo_b.size()];

   for(int elec1=0;elec1<int(mo_a.size());elec1++)
   {
      for(int elec2=0;elec2<int(mo_b.size());elec2++)
      {
         if(spin_a.at(elec1) == spin_b.at(elec2))
         matrix[elec1*mo_b.size()+elec2]=MO_S[mo_a.at(elec1)*n_occ+mo_b.at(elec2)];
      }
   }

   sum=determinant(matrix,mo_b.size());

   delete [] matrix;

   return sum;

}
