// This is a compilation of functions extracted from S2kit 1.0
// done by Henrique S. Xavier (hsxavier@if.usp.br) on Nov-2014.
// Only functions needed to perform Forward and Inverse Discrete 
// Legendre Transforms (DLT) were selected.


/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


#include <math.h>
#include <string.h>   /* to declare memcpy */
#include <stdlib.h>
#include <stdio.h>
#ifndef PI
#define PI 3.14159265358979
#endif
#include "interpol.h"
#include "Utilities.hpp"
#include "s2kit10_naive.hpp"


/************************************************************************/
/* returns an array of the angular arguments of n Chebyshev nodes */
/* eval_pts points to a double array of length n */

void ArcCosEvalPts(int n,
		   double *eval_pts)
{
    int i;
    double twoN;

    twoN = (double) (2 * n);

   for (i=0; i<n; i++)
     eval_pts[i] = (( 2.0*((double)i)+1.0 ) * PI) / twoN;

}
/************************************************************************/
/* returns an array of n Chebyshev nodes */

void EvalPts( int n,
	      double *eval_pts)
{
    int i;
    double twoN;

    twoN = (double) (2*n);

   for (i=0; i<n; i++)
     eval_pts[i] = cos((( 2.0*((double)i)+1.0 ) * PI) / twoN);

}


/************************************************************************/
/* vector arithmetic operations */
/************************************************************************/
/* does result = data1 + data2 */
/* result and data are vectors of length n */

void vec_add(double *data1,
	     double *data2,
	     double *result,
	     int n)
{
  int k;


  for (k = 0; k < n % 4; ++k)
    result[k] = data1[k] + data2[k];

  for ( ; k < n ; k += 4)
    {
      result[k] = data1[k] + data2[k];
      result[k + 1] = data1[k + 1] + data2[k + 1];
      result[k + 2] = data1[k + 2] + data2[k + 2];
      result[k + 3] = data1[k + 3] + data2[k + 3];
    }
}
/************************************************************************/
/************************************************************************/
/*
   vec_mul(scalar,data1,result,n) multiplies the vector 'data1' by
   'scalar' and returns in result 
*/
void vec_mul(double scalar,
	     double *data1,
	     double *result,
	     int n)
{
   int k;


   for( k = 0; k < n % 4; ++k)
     result[k] = scalar * data1[k];

   for( ; k < n; k +=4)
     {
       result[k] = scalar * data1[k];
       result[k + 1] = scalar * data1[k + 1];
       result[k + 2] = scalar * data1[k + 2];
       result[k + 3] = scalar * data1[k + 3];
     }

}


/************************************************************************/
/* point-by-point multiplication of vectors */

void vec_pt_mul(double *data1,
		double *data2,
		double *result,
		int n)
{
   int k;

  
  for(k = 0; k < n % 4; ++k)
    result[k] = data1[k] * data2[k];
  
  for( ; k < n; k +=4)
    {
      result[k] = data1[k] * data2[k];
      result[k + 1] = data1[k + 1] * data2[k + 1];
      result[k + 2] = data1[k + 2] * data2[k + 2];
      result[k + 3] = data1[k + 3] * data2[k + 3];
    }
 
}


/************************************************************************/
/* Recurrence coefficients */
/************************************************************************/
/* Recurrence coefficents for L2-normed associated Legendre
   recurrence.  When using these coeffs, make sure that
   inital Pmm function is also L2-normed */
/* l represents degree, m is the order */

double L2_an(int m,
	     int l)
{
  return (sqrt((((double) (2*l+3))/((double) (2*l+1))) *
	       (((double) (l-m+1))/((double) (l+m+1)))) *
	  (((double) (2*l+1))/((double) (l-m+1))));

}

/* note - if input l is zero, need to return 0 */
double L2_cn(int m,
	     int l) 
{
  if (l != 0) {
    return (-1.0 *
	  sqrt((((double) (2*l+3))/((double) (2*l-1))) *
	       (((double) (l-m+1))/((double) (l+m+1))) *
	       (((double) (l-m))/((double) (l+m)))) *
	  (((double) (l+m))/((double) (l-m+1))));
  }
  else
    return 0.0;

}


/************************************************************************/
/* L2 normed Pmm.  Expects input to be the order m, an array of
 evaluation points arguments of length n, and a result vector of length n */
/* The norming constant can be found in Sean's PhD thesis */
/* This has been tested and stably computes Pmm functions thru bw=512 */

void Pmm_L2( int m,
	     double *eval_pts,
	     int n,
	     double *result)
{
  int i;
  double md, id, mcons;

  id = (double) 0.0;
  md = (double) m;
  mcons = sqrt(md + 0.5);

  for (i=0; i<m; i++) {
    mcons *= sqrt((md-(id/2.0))/(md-id));
    id += 1.0;
  }
  if (m != 0 )
    mcons *= pow(2.0,-md/2.0);
  if ((m % 2) != 0) mcons *= -1.0;

  for (i=0; i<n; i++) 
    result[i] = mcons * pow(sin(eval_pts[i]),((double) m));

}


/************************************************************************/
/* generate all of the Pmls for a specified value of m.  

   storeplm points to a double array of size 2 * bw * (bw - m);

   Workspace needs to be
   16 * bw 

   P(m,l,j) respresents the associated Legendre function P_l^m
   evaluated at the j-th Chebyshev point (for the bandwidth bw)
   Cos((2 * j + 1) * PI / (2 * bw)).

   The array is placed in storeplm as follows:

   P(m,m,0)    P(m,m,1)  ... P(m,m,2*bw-1)
   P(m,m+1,0)  P(m,m+1,1)... P(m,m+1,2*bw-1)
   P(m,m+2,0)  P(m,m+2,1)... P(m,m+2,2*bw-1)
   ...
   P(m,bw-1,0)   P(m,bw-1,1) ... P(m,bw-1,2*bw-1)

   This array will eventually be used by the naive transform algorithm.
   This function will precompute the arrays necessary for the algorithm.
*/
void PmlTableGen(int bw,
		 int m,
		 double *storeplm,
		 double *workspace)
{
  double *prev, *prevprev;
  double *temp1, *temp2, *temp3, *temp4;
  double *x_i, *eval_args;
  int i;
  
  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);
  

  /* get the evaluation nodes */
  EvalPts(2*bw,x_i);
  ArcCosEvalPts(2*bw,eval_args);
  
  /* set initial values of first two Pmls */
  for (i=0; i<2*bw; i++) 
    prevprev[i] = 0.0;
  if (m == 0)
    for (i=0; i<2*bw; i++)
      prev[i] = 0.707106781186547 ;
  else 
    Pmm_L2(m, eval_args, 2*bw, prev);

  memcpy(storeplm, prev, sizeof(double) * 2 * bw);

  for(i = 0; i < bw - m - 1; i++)
    {
      vec_mul(L2_cn(m,m+i),prevprev,temp1,2*bw);
      vec_pt_mul(prev, x_i, temp2, 2*bw);
      vec_mul(L2_an(m,m+i), temp2, temp3, 2*bw);
      vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
      
      storeplm += (2 * bw);
      memcpy(storeplm, temp4, sizeof(double) * 2 * bw);
      memcpy(prevprev, prev, sizeof(double) * 2 * bw);
      memcpy(prev, temp4, sizeof(double) * 2 * bw);
    }
}


/************************************************************************/
/* This is the procedure that synthesizes a function from a list
   of coefficients of a Legendre series. I.e. this is the INVERSE
   discrete Legendre transform.

   Function is synthesized at the (2*bw) Chebyshev nodes. Associated
   Legendre functions are assumed to be precomputed.
   
   bw - bandwidth

   m - order

   coeffs - a pointer to double array of size (bw-m).  First coefficient is
            coefficient for Pmm
   result - a pointer to double array of size (2*bw); at the conclusion
            of the routine, this array will contain the
            synthesized function

   plmtable - a pointer to a double array of size (2*bw*(bw-m));
	      contains the PRECOMPUTED plms, i.e. associated Legendre
	      functions. E.g. Should be generated by a call to

	      PmlTableGen(),

	      (see pmls.c)


	      NOTE that these Legendres are normalized with norm
	      equal to 1 !!!
*/

void Naive_SynthesizeX(double *coeffs,
		       int bw,
		       int m,
		       double *result,
		       double *plmtable)
{
  int i, j;
  double tmpcoef;

  /* make sure result is zeroed out */
  memset( result, 0, sizeof(double) * 2 * bw );

  for ( i = 0 ; i < bw - m ; i ++ )
    {
      tmpcoef = coeffs[i] ;
      if ( tmpcoef != 0.0 )
	for (j=0; j<(2*bw); j++)
	  result[j] += (tmpcoef * plmtable[j]);
      plmtable += (2 * bw ) ;
    }
}

/*
  makeweights: given a bandwidth bw, make weights for
  both even *and* odd order Legendre transforms.

  bw -> bandwidth of transform
  weights -> pointer to double array of length 4*bw; this
             array will contain the even and odd weights;
	     even weights start at weight[0], and odd weights
	     start at weights[2*bw]

*/

void makeweights( int bw,
		  double *weights )
{
  int j, k ;
  double fudge ;
  double tmpsum ;

  fudge = M_PI/((double)(4*bw)) ;
  

  for ( j = 0 ; j < 2*bw ; j ++ )
    {
      tmpsum = 0.0 ;
      for ( k = 0 ; k < bw ; k ++ )
	tmpsum += 1./((double)(2*k+1)) *
	  sin((double)((2*j+1)*(2*k+1))*fudge);
      tmpsum *= sin((double)(2*j+1)*fudge);
      tmpsum *= 2./((double) bw) ;
      
      weights[j] = tmpsum ;
      weights[j + 2*bw] = tmpsum * sin((double)(2*j+1)*fudge);
    }

}


/************************************************************************/
/*

   Naive_AnalysisX: computing the discrete Legendre transform of
                    a function via summing naively. I.e. This is
		    the FORWARD discrete Legendre transform.

   bw - bandwidth
   m - order

   data - a pointer to double array of size (2*bw) containing
          the sample points

   result - a pointer to double array of size (bw-m) which, at the
            conclusion of the routine, will contains the coefficients

   plmtable - a pointer to a double array of size (2*bw*(bw-m));
	      contains the PRECOMPUTED plms, i.e. associated Legendre
	      functions.  E.g. Should be generated by a call to

	      PmlTableGen()

	      (see pmls.c)

	      NOTE that these Legendres are normalized with norm
	      equal to 1 !!!

   workspace - array of size 2 * bw;


*/


void Naive_AnalysisX(double *data,
		     int bw,
		     int m,
		     double *weights,
		     double *result,
		     double *plmtable,
		     double *workspace)
{
  int i, j;
  double result0, result1, result2, result3;
  register double *wdata;

  wdata = workspace;

  /* make sure result is zeroed out */
  memset( result, 0, sizeof(double) * (bw - m) );

  /* apply quadrature weights */
  /*
    I only have to differentiate between even and odd
    weights when doing something like seminaive, something
    which involves the dct. In this naive case, the parity of
    the order of the transform doesn't matter because I'm not
    dividing by sin(x) when precomputing the Legendres (because
    I'm not taking their dct). The plain ol' weights are just
    fine. 
  */
  
  for(i = 0; i < 2 * bw; i++)
    wdata[i] = data[i] * weights[i];

  /* unrolling seems to work */
  if ( 1 )
    {
      for (i = 0; i < bw - m; i++)
	{
	  result0 = 0.0; result1 = 0.0;
	  result2 = 0.0; result3 = 0.0;
	  
	  for(j = 0; j < (2 * bw) % 4; ++j)
	    result0 += wdata[j] * plmtable[j];

	  for( ; j < (2 * bw); j += 4)
	    {
	      result0 += wdata[j] * plmtable[j];
	      result1 += wdata[j + 1] * plmtable[j + 1];
	      result2 += wdata[j + 2] * plmtable[j + 2];
	      result3 += wdata[j + 3] * plmtable[j + 3];
	    }
	  result[i] = result0 + result1 + result2 + result3;

	  plmtable += (2 * bw);
	}
    }
  else
    {
      for (i = 0; i < bw - m; i++)
	{
	  result0 = 0.0 ; 	  
	  for(j = 0; j < (2 * bw) ; j++)
	    result0 += wdata[j] * plmtable[j];
	  result[i] = result0 ;

	  plmtable += (2 * bw);
	}
    }
}


/******************************************************************/
/****** Functions written by Henrique S. Xavier to ease     *******/
/****** use of s2kit10.                                     *******/
/******************************************************************/


/*** Exponential suppression of C_l to avoid oscillations in Xi ***/
double suppress(double l, double lsup, double supindex) {
  return exp( -1.0*pow(l/lsup, supindex) );
}


/*** Generate C_l coefficients for Discrete Legendre Transform implemented by s2kit10 ***/
/*** Also interpolates the input Cls to get a C_l for each l                          ***/
double *GetCl4DLT(double *Clin, double *ll, int Clinsize, double lsup, double supindex, int lmax) {
  const double sqr2over4pi = 0.1125395395196383;
  int l;
  double *Clout;

  Clout = vector<double>(0, lmax);
  for(l=1; l<=lmax; l++) 
    Clout[l] = sqr2over4pi*sqrt((double)(2*l+1))*Interpol(ll, Clinsize, Clin, (double)l)*suppress((double)l,lsup,supindex);
  Clout[0] = 0.0; // Set C_0=0.

  return Clout;
}


/*** Interpolates the input Cls to get a C_l for each l, from 0 to lmax ***/
void GetAllLs(double *ll, double *Clin, int Clinsize, double *Clout, int lmax, int extrapol) {
  int l;

                                Clout[0] = 0.0; // Set monopole to zero.
  if (ll[0]>1.0 && extrapol==0) Clout[1] = 0.0; // If not present, set dipole to zero too.
  else                          Clout[1] = Interpol(ll, Clinsize, Clin, 1.0);
  for(l=2; l<=lmax; l++)        Clout[l] = Interpol(ll, Clinsize, Clin, (double)l);
}


/*** Transform C_l coefficients for Discrete Legendre Transform implemented by s2kit10 ***/
void ModCl4DLT(double *Cl, int lmax, double lsup, double supindex) {
  const double sqr2over4pi = 0.1125395395196383;
  int l;
  
  // Suppression needed to avoid oscillations in Xi(theta).
  if (lsup>0 && supindex>0) {
    for(l=0; l<=lmax; l++) 
      Cl[l] = sqr2over4pi*sqrt((double)(2*l+1))*Cl[l]*suppress((double)l,lsup,supindex);
  }
  // Suppression can be turned off by seeting lsup or supindex to <0.
  else {
    for(l=0; l<=lmax; l++) 
      Cl[l] = sqr2over4pi*sqrt((double)(2*l+1))*Cl[l];
  }
}


/*** Erase suppression generated by function 'suppress' ***/
double unsuppress(double l, double lsup, double supindex) {
  return exp( pow(l/lsup, supindex) );
}


/*** Apply factors to Cl obtained from s2kit10 routines ***/
/*** May unsuppress suppression, but only if Cls are for the same distribution as before ***/
void ApplyClFactors(double *Cl, int ClLength, double lsup/*=-1*/, double supindex/*=0*/) {
  const double twopisqr2 = 8.885765876316732; // 2*PI*sqrt(2).
  int l;
  
  if (lsup<0 || supindex<0) for (l=0; l<ClLength; l++) Cl[l] = twopisqr2*Cl[l]/sqrt((double)(2*l+1));
  else for (l=0; l<ClLength; l++) Cl[l] = twopisqr2*Cl[l]/sqrt((double)(2*l+1))*unsuppress((double)l, lsup, supindex);
}
