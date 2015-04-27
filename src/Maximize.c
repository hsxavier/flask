#include "interpol.h"

/*** Code based on functions from Numerical Recipes mnbrak and golden.             ***/
/*** Modified by Henrique S. Xavier in 24/apr/2015 to be used with interpolations. ***/
/*** Modified also to look for maximums instead of minimums.                       ***/


// From nrutil.h:
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
		      (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
// Definitions from Numerical Recipes:
#include <math.h>
#define NRANSI
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


// Interpolation specifically for use with functions below.
// Returns -f(x) if xmin<x<xmax or infinite otherwise. Made for searches of maximum within xmin and xmax.
double BoundInterpol(double *Xarray,int NX,double *Yarray,double x, double xmin, double xmax) {
  if(x<xmin || x>xmax) return 1.0/0.0;
  else return -1.0*Interpol(Xarray,NX,Yarray, x);
}

/*** Functions from Numerical Recipes, modified ***/

// Function to bracket a minimum, stay inside boundaries:
void mnbrakInterp(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, // << mnbrak function parameters.
		  double *Xarray,int NX,double *Yarray) // << Interpolation parameters.
{
  double ulim,u,r,q,fu,dum, xmin, xmax;

  xmin = *ax;  xmax = *bx;
  *fa=BoundInterpol(Xarray,NX,Yarray, *ax, xmin, xmax);
  *fb=BoundInterpol(Xarray,NX,Yarray, *bx, xmin, xmax);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=BoundInterpol(Xarray,NX,Yarray, *cx, xmin, xmax);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=BoundInterpol(Xarray,NX,Yarray, u, xmin, xmax);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */


#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);


// Find the minimum already bracketed, stay inside boudaries:
double goldenInterp(double ax, double bx, double cx, double tol, double *xmin, // << golden function parameters.
		    double *Xarray,int NX,double *Yarray, double xmin0, double xmax0) // << Interpolation parameters.
{
  double f1,f2,x0,x1,x2,x3;

  x0=ax;
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax)) {
    x1=bx;
    x2=bx+C*(cx-bx);
  } else {
    x2=bx;
    x1=bx-C*(bx-ax);
  }
  f1=BoundInterpol(Xarray,NX,Yarray, x1, xmin0, xmax0);
  f2=BoundInterpol(Xarray,NX,Yarray, x2, xmin0, xmax0);
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) {
      SHFT3(x0,x1,x2,R*x1+C*x3)
	SHFT2(f1,f2,BoundInterpol(Xarray,NX,Yarray, x2, xmin0, xmax0))
	} else {
      SHFT3(x3,x2,x1,R*x2+C*x0)
	SHFT2(f2,f1,BoundInterpol(Xarray,NX,Yarray, x1, xmin0, xmax0))
	}
  }
  if (f1 < f2) {
    *xmin=x1;
    return f1;
  } else {
    *xmin=x2;
    return f2;
  }
}
#undef C
#undef R
#undef SHFT2
#undef SHFT3
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */


// Returns the maximum of an interpolated function within the intervals xmin and xmax:
double MaxInterp(double xmin, double xmax, double tol, double *Xarray,int NX,double *Yarray) {
  double ax, bx, cx, fa, fb, fc, maximum, xmaximum;
  
  ax = xmin; bx = xmax;
  mnbrakInterp(&ax, &bx, &cx, &fa, &fb, &fc, Xarray, NX, Yarray);
  maximum = goldenInterp(ax, bx, cx, tol, &xmaximum, Xarray, NX, Yarray, xmin, xmax);
  return -1.0*maximum;

}
