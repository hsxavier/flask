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


/*** Compute the integral of a tabulated function by summing trapeziums ***/
/*** I.e., this uses linear interpolation between points.               ***/
double DiscreteIntegral(double *x, double *f, int n, double xmin, double xmax) {
  int i, imin, imax;
  double integral=0.0;
  
  // Guard against integrating beyond tabulated values:
  if (xmin<x[0])   { printf("DiscreteIntegral ERROR: xmin<x[0]\n"  ); exit(-1); }
  if (xmax>x[n-1]) { printf("DiscreteIntegral ERROR: xmax>x[n-1]\n"); exit(-1); }
  
  // Find points that are wrapped by xmin and xmax:
  bisect_interpol(x, n, xmin, &imin); imin=imin+1;
  bisect_interpol(x, n, xmax, &imax);
  
  // Integrate over wrapped points:
  for (i=imin; i<imax; i++) integral += (f[i]+f[i+1]) * (x[i+1]-x[i]) / 2.0;
  // Add contribution of borders (interpolated points):
  integral += (Interpol(x,n,f,xmin)+f[imin]) * (x[imin]-xmin) / 2.0;
  integral += (Interpol(x,n,f,xmax)+f[imax]) * (xmax-x[imax]) / 2.0;

  return integral;
}
