#include "Spline.hpp"
#include "Utilities.hpp"
#include <cstddef> // For NULL pointer.

// Originally the Numerical recipes functions were written for vectors x[1..n] (i.e. first element at i=1).
// The constant below controls the offset of the vector consistently across the functions in this file.
// Now vectors are x[offset..n+offset-1] 
// The comments before each Num-Rec function were not updated.
const int splOffset = 0;
const double deriv0 = 1e40; // If this is >1e30, the second derivatives at edges are set to zero.



/*************************************/
/***         Spline class          ***/
/*************************************/


/*** Empty constructor ***/
Spline::Spline() {
  clear();
}


/*** Construct and initialize 1D spline ***/
Spline::Spline(double *dom, double *f, int n) {
  clear();
  init(dom, f, n);
}


/*** Construct and initialize 2D spline ***/
Spline::Spline(double *dom1, double *dom2, double **f, int n1, int n2) {
  clear();
  init(dom1, dom2, f, n1, n2);
}


/*** Initialize 1D spline ***/
void Spline::init(double *x, double *y, int n) {
  int i;
  
  // Clear previous settings if existent:
  if (dim!=0) dealloc();
  // Set variables and allocate memory:
  dim     = 1;
  N1      = n;
  N2      = 0;
  x1      = vector<double>(splOffset, N1+splOffset-1);
  x2      = NULL;
  f1d     = vector<double>(splOffset, N1+splOffset-1);
  f2d     = NULL;
  d2fdx2  = vector<double>(splOffset, N1+splOffset-1);
  d2fdxdy = NULL;
  // Copy function sampled point to internal arrays:
  for (i=splOffset; i<N1+splOffset; i++) {
    x1 [i] = x[i];
    f1d[i] = y[i];
  }
  // Compute second derivatives that are needed for Spline:
  spline(x1, f1d, N1, deriv0, deriv0, d2fdx2);
}


/*** Initialize 2D spline ***/
void Spline::init(double *dom1, double *dom2, double **f, int n1, int n2) {
  int i, j;
  
  // Clear previous settings if existent:
  if (dim!=0) dealloc();
  // Set variables and allocate memory:
  dim     = 2;
  N1      = n1;
  N2      = n2;
  x1      = vector<double>(splOffset, N1+splOffset-1);
  x2      = vector<double>(splOffset, N2+splOffset-1);
  f1d     = NULL;
  f2d     = matrix<double>(splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
  d2fdx2  = NULL;
  d2fdxdy = matrix<double>(splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
  // Copy function sampled point to internal arrays:
  for (i=splOffset; i<N1+splOffset; i++) x1[i] = dom1[i];
  for (j=splOffset; j<N2+splOffset; j++) x2[j] = dom2[j];
  for (i=splOffset; i<N1+splOffset; i++) 
    for (j=splOffset; j<N2+splOffset; j++){
      f2d[i][j] = f[i][j];
  }
  // Compute second derivatives that are needed for Spline:
  splie2(x1, x2, f2d, N1, N2, d2fdxdy);
}


/*** Loading specifically for the flask and Dens2KappaCls program ***/
// Will use the Cov[i][j][l] tensor, you have to specify the field index and 
// the multipole l, and the code loops over redshifts.
void Spline::init(const FZdatabase & fieldlist, double ***Cov, int field, int lpos) {
  int i, j, z1, z2;
  
  // Clear previous settings if existent:
  if (dim!=0) dealloc();
  // Set variables and allocate memory:
  dim     = 2;
  N1      = fieldlist.Nz4f(field);
  N2      = fieldlist.Nz4f(field);
  x1      = vector<double>(splOffset, N1+splOffset-1);
  x2      = vector<double>(splOffset, N2+splOffset-1);
  f1d     = NULL;
  f2d     = matrix<double>(splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
  d2fdx2  = NULL;
  d2fdxdy = matrix<double>(splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
  // Copy function sampled point to internal arrays:
  for (z1=0; z1<N1; z1++) {
    fieldlist.fFixedIndex(field, z1, &i);
    // Will use the mean redshift in the bin as point where Cov was sampled:
    x1[splOffset+z1] = (fieldlist.zmin(i)+fieldlist.zmax(i))/2.0;
    x2[splOffset+z1] = x1[splOffset+z1];
  }
  for (z1=splOffset; z1<N1+splOffset; z1++) {
    fieldlist.fFixedIndex(field, z1, &i);
    for (z2=splOffset; z2<N2+splOffset; z2++) { 
      fieldlist.fFixedIndex(field, z2, &j);
      f2d[z1][z2] = Cov[i][j][lpos];
    }
  }
  // Compute second derivatives that are needed for Spline:
  splie2(x1, x2, f2d, N1, N2, d2fdxdy);
}


/*** Returns the 1D cubic spline interpolated value ***/
double Spline::operator()(double xpt) const {
  double result;
  if (dim!=1) error("Spline.operator(): wrong dimension (expected 1) or not initialized.");
  splint(x1, f1d, d2fdx2, N1, xpt, &result);
  return result;
}


/*** Returns the 2D cubic spline interpolated value ***/
double Spline::operator()(double x1pt, double x2pt) const {
  double result;
  if (dim!=2) error("Spline.operator(): wrong dimension (expected 2) or not initialized.");
  splin2(x1, x2, f2d, d2fdxdy, N1, N2, x1pt, x2pt, &result);
  return result;
}


/*** Destructor, deallocate internal memory ***/
Spline::~Spline() {
  dealloc();
}


/*** Deallocate internal memory ***/
void Spline::dealloc() {
  if (dim!=0) {
    if      (dim==1) {
      free_vector(x1,      splOffset, N1+splOffset-1);
      free_vector(f1d,     splOffset, N1+splOffset-1);
      free_vector(d2fdx2,  splOffset, N1+splOffset-1); 
    }
    else if (dim==2) {
      free_vector(x1,      splOffset, N1+splOffset-1);
      free_vector(x2,      splOffset, N2+splOffset-1);
      free_matrix(f2d,     splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
      free_matrix(d2fdxdy, splOffset, N1+splOffset-1, splOffset, N2+splOffset-1);
    }
    else warning("Spline.dealloc: unknown dimension, not deallocating.");
  }
}


/*** Clear internal variables ***/
void Spline::clear() {
  dim     = 0;
  N1      = 0;
  N2      = 0;
  x1      = NULL;
  x2      = NULL;
  f1d     = NULL;
  f2d     = NULL;
  d2fdx2  = NULL;
  d2fdxdy = NULL;
}


/*** Fully clear object ***/
void Spline::reset() {
  dealloc();
  clear();
}



/*************************************/
/*** Standalone backbone functions ***/
/*************************************/

#define NRANSI

/*
Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., y i = f ( x i ), with
x 1 < x 2 < . . . < x N , and given values yp1 and ypn for the first derivative of the interpolating
function at points 1 and n , respectively, this routine returns an array y2[1..n] that contains
the second derivatives of the interpolating function at the tabulated points x i . If yp1 and/or
ypn are equal to 1 × 10 30 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary.
*/

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]) {
  int i,k;
  double p,qn,sig,un,*u;

  u=vector<double>(splOffset,n+splOffset-2);
  if (yp1 > 0.99e30)
    y2[splOffset]=u[splOffset]=0.0;
  else {
    y2[splOffset] = -0.5;
    u[splOffset]=(3.0/(x[splOffset+1]-x[splOffset]))*((y[splOffset+1]-y[splOffset])/(x[splOffset+1]-x[splOffset])-yp1);
  }
  for (i=splOffset+1;i<n+splOffset-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n+splOffset-1]-x[n+splOffset-2]))*(ypn-(y[n+splOffset-1]-y[n+splOffset-2])/(x[n+splOffset-1]-x[n+splOffset-2]));
  }
  y2[n+splOffset-1]=(un-qn*u[n+splOffset-2])/(qn*y2[n+splOffset-2]+1.0);
  for (k=n+splOffset-2;k>=splOffset;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,splOffset,n+splOffset-2);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */


/*
Given the arrays xa[1..n] and ya[1..n] , which tabulate a function (with the xa i ’s in order),
and given the array y2a[1..n] , which is the output from spline above, and given a value of
x , this routine returns a cubic-spline interpolated value y .
*/
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y) {
  int klo,khi,k;
  double h,b,a;

  klo=splOffset;
  khi=n+splOffset-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) error("splint: Bad xa input to routine splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */


/*
Given an m by n tabulated function ya[1..m][1..n] , and tabulated independent variables
x2a[1..n] , this routine constructs one-dimensional natural cubic splines of the rows of ya
and returns the second-derivatives in the array y2a[1..m][1..n] . (The array x1a[1..m] is
included in the argument list merely for consistency with routine splin2 .)
*/
void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a) {
  int j;
  for (j=splOffset;j<m+splOffset;j++) spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */



/*
Given x1a , x2a , ya , m , n as described in splie2 and y2a as produced by that routine; and
given a desired interpolating point x1 , x2 ; this routine returns an interpolated function value y
by bicubic spline interpolation.
*/

#define NRANSI
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
	    double x1, double x2, double *y) {
  int j;
  double *ytmp,*yytmp;
  
  ytmp=vector<double>(splOffset,m+splOffset-1);
  yytmp=vector<double>(splOffset,m+splOffset-1);
  for (j=splOffset;j<m+splOffset;j++) splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  splint(x1a,yytmp,ytmp,m,x1,y);
  free_vector(yytmp,splOffset,m+splOffset-1);
  free_vector(ytmp,splOffset,m+splOffset-1);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */

