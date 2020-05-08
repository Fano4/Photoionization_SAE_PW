#include <cmath>

double phi(int m,double x)
{
   if(m<0)
      return sin(abs(m)*x);
   else
      return cos(m*x);

}
