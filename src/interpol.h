#ifndef INTERPOL_H
#define INTERPOL_H 1

#include <stdio.h>
#include <stdlib.h>


/*** Part of interpolation function below ***/
/*** From Benjamin Joachimi ***/
void bisect_interpol(double *a,int dim,double x,int *erg);

/*** Interpolation of function with domain ax and values z, NX points, at x. ***/
/*** Originally 'interpol_linear_extra_1D' from Benjamin Joachimi's code. ***/
double Interpol(double *ax,int NX,double *z,double x);

/*** Compute the integral of a tabulated function by summing trapeziums ***/
/*** I.e., this uses linear interpolation between points.               ***/
double DiscreteIntegral(double *x, double *f, int n, double xmin, double xmax);

#endif
