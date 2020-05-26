#include <cstdlib>
#include <iostream>
#include <vector>

void sym_to_nosym_mat_trans(int n_sym,std::vector<int> dim_1_sym,std::vector<int> dim_2_sym,std::vector<double> mat_sym_rep,int* dim1,int* dim2,std::vector<double>* mat_nosym_rep)
{
   //Transforms a block matrix written in symmetry-consistent representation to a single block matrix
   //The mat_sym_rep matrix has a dimensionality sum_i ( dim1[i]*dim2[i] )
   //The mat_nosym_rep matrix has a dimensionality (sum_i dim1[i])*(sum_i dim2[i]) 
   

   //compute total matrix size
   *dim1=0;
   *dim2=0;
   for(int i=0;i!=n_sym;i++)
   {
      *dim1+=dim_1_sym.at(i);
      *dim2+=dim_2_sym.at(i);
   }

   //initialize the mat_nosym_rep vector
   mat_nosym_rep->clear();
   for(int i=0;i!=*dim1**dim2;i++)
      mat_nosym_rep->push_back(0);

   //Tranform the representation
   int index(0);
   int index1(0);
   int index2(0);
   int mem2(0);
   for(int i=0;i!=n_sym;i++)
   {
      for(int c=0;c!=dim_1_sym[i];c++)
      {
         index2=mem2;
         for(int s=0;s!=dim_2_sym[i];s++)
         {
            mat_nosym_rep->at(index1**dim2+index2)=mat_sym_rep.at(index);
            index++;
            index2++;
         }
         index1++;
      }
      mem2=index2;
   }

}
