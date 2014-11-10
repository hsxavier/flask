#include <stdio.h>
#include <stdlib.h>
#include "interpol.h"

/*** Part of interpolation function below ***/
/*** From Benjamin Joachimi ***/
void bisect_interpol(double *a,int dim,double x,int *erg)
{
  int l=0,u=dim-1,m;

  if ((x<a[l]) || (x>a[u])) {
    printf("Error: index out of range: %g not in %g-%g\n",x,a[l],a[u]);
    exit(-1);
  }
  while (u-l>1) {
    m=(u+l)/2;
    if (x>a[m]) l=m;
    else u=m;
  }
  *erg=l;
}


/*** Interpolation of function with domain ax and values z, NX points, at x. ***/
/*** Originally 'interpol_linear_extra_1D' from Benjamin Joachimi's code. ***/
double Interpol(double *ax,int NX,double *z,double x) {

  if (x<=ax[0]) return((z[1]-z[0])/(ax[1]-ax[0])*(x-ax[0])+z[0]);
  else if (x>=ax[NX-1]) return((z[NX-1]-z[NX-2])/(ax[NX-1]-ax[NX-2])*(x-ax[NX-1])+z[NX-1]);
  else {
    int ix;
    bisect_interpol(ax,NX,x,&ix);
    double t=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    return((1-t)*z[ix]+t*z[ix+1]);
  }
}
