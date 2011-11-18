/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <math.h>
#include <iostream>
using namespace std;

/* ----------------------------------------------------------------------
   Evaluate log(Gamma(x))
   From Tom Minka (Microsoft Research Cambridge)
---------------------------------------------------------------------- */

double gammaln(double x)
{
  const static double M_lnSqrt2PI=0.91893853320467274178;
  static double gamma_series[] = {
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  };
  int i;
  double denom, x1, series;
  if(x < 0) return 1e308;
  denom = x+1;
  x1 = x + 5.5;
  series = 1.000000000190015;
  for(i = 0; i < 6; i++) {
    series += gamma_series[i] / denom;
    denom += 1.0;
  }
  return M_lnSqrt2PI + (x+0.5)*log(x1) - x1 + log(series/x);
}

/* ----------------------------------------------------------------------
   This is a fast version of the log(Gamma(x)) function.
   
   It precomputes a lookup table for the interval [LB,UB2] on the first
   call, and uses quadratic interpolation to approximate gammaln on this
   interval.

   It is typically accurate to 5 d.p. and stores 48 KB on the stack.
---------------------------------------------------------------------- */

double fast_gammaln(double x)
{
  const static double LB=0.01;
  const static double UB=10;
  const static double LB2=UB;
  const static double UB2=50;
  const static int GRANULARITY=4096;
  const static int GRANULARITY2=2048;
  static double lookup[GRANULARITY];
  static double lookup2[GRANULARITY2];
  static bool lookup_init=false;
  double X,x1,x2,x3,y1,y2,y3,A,B,C;
  int X1,X2,X3;
  if(x<LB || x>UB2) return gammaln(x);
  if(!lookup_init)
    {
      double inc=(UB-LB)/(double)GRANULARITY;
      int id;
      X=LB;
      for(id=0; id<GRANULARITY; id++)
	{
	  lookup[id]=gammaln(X);
	  X+=inc;
	}
      inc=(UB2-LB2)/(double)GRANULARITY2;
      X=LB2;
      for(id=0; id<GRANULARITY2; id++)
	{
	  lookup2[id]=gammaln(X);
	  X+=inc;
	}
      lookup_init=true;
    }
  if(x>UB)
    {
      X=(x-LB2)*(double)GRANULARITY2/(UB2-LB2);
      X1=(int)X;
      X2=X1+1;
      X3=X2+1;
      while(X3 >= GRANULARITY2)
	{
	  X1--; X2--; X3--;
	}
      x1=X1; x2=X2; x3=X3;
      y1=lookup2[X1];
      y2=lookup2[X2];
      y3=lookup2[X3];
    }else
    {
      X=(x-LB)*(double)GRANULARITY/(UB-LB);
      X1=(int)X;
      X2=X1+1;
      X3=X2+1;
      while(X3 >= GRANULARITY)
	{
	  X1--; X2--; X3--;
	}
      x1=X1; x2=X2; x3=X3;
      y1=lookup[X1];
      y2=lookup[X2];
      y3=lookup[X3];
    }

  A=((y3-y1)+((x3-x1)/(x2-x1))*(y2-y1))/
    ((x3*x3-x1*x1)+((x3-x1)/(x2-x1))*(x2*x2-x1*x1));
  B=((y2-y1)-A*(x2*x2-x1*x1))/(x2-x1);
  C=y1-A*x1*x1-B*x1;
  
  return A*X*X + B*X + C;
}
