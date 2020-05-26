#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include "utilities.h"
#include "mathfunctions.h"

void slater_ovlp(int n_elec,int n_csf_a,int n_csf_b,std::vector<int> csf_mo_a,std::vector<int> csf_mo_b,std::vector<int> csf_spin_a,std::vector<int> csf_spin_b,int n_occ_a,int n_occ_b,std::vector<double> MO_S,std::vector<double>* CSF_S)
{
   //This function computes the overlap between two slater determinants without computing a determinant. 
   //This version requires the two determinants to have the same size
   
   using namespace std;

   double sum(0);

   CSF_S->clear();
   double* matrix=new double[n_elec*n_elec];
   
   for(int csfa=0;csfa<n_csf_a;csfa++)
   {
      for(int csfb=0;csfb<n_csf_b;csfb++)
      {

         for(int elec1=0;elec1<n_elec;elec1++)
         {
            for(int elec2=0;elec2<n_elec;elec2++)
            {
               if(csf_spin_a.at(csfa*n_elec+elec1) == csf_spin_b.at(csfb*n_elec+elec2))
                  matrix[elec1*n_elec+elec2]=MO_S.at(csf_mo_a.at(csfa*n_elec+elec1)*n_occ_b+csf_mo_b.at(csfb*n_elec+elec2));
            }
         }
         CSF_S->push_back(determinant(matrix,n_elec));
      }
   }
   delete [] matrix;

}
