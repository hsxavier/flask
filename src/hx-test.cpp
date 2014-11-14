#include <stdio.h>
//#include "naive_synthesis.h"
//#include "pmls.h"
//#include "primitive.h"
//#include "naive_synthesis.c"
//#include "pmls.c"
//#include "primitive.c"
#include "s2kit10_naive.h"
#include "../src/Utilities.hpp"
#include "../src/interpol.h"           // Interpolation. 
#include <math.h>

int main () {
  double *ll, *Cov, *wrapper[2], *workspace, *plm, *coeffs, *xi, *theta, *weights, *newcoeffs; 
  int bw=10000, ClLength=200, maxL=bw, l;
  double factor;

  // Load C_l:
  printf("Load input Cl.\n");
  ll  = vector<double>(0, ClLength-1);
  Cov = vector<double>(0, ClLength-1);
  wrapper[0]=ll;
  wrapper[1]=Cov;
  ImportVecs(wrapper, ClLength, 2, "../data/ps_test_grf.dat");
  coeffs = GetCl4DLT(Cov, ll, ClLength, 8000, 6.0, maxL);
  //for(l=0; l<maxL; l++) coeffs[l] = sqrt((double)(2*l+1))*Interpol(ll, ClLength, Cov, (double)l);
  //coeffs[0]=0; // Set C_0 = 0
  
  // Prepare for trafo
  printf("Generate plms.\n");
  workspace = vector<double>(0, 16*bw-1);
  plm = vector<double>(0, 2*bw*bw-1);
  PmlTableGen(bw, 0, plm, workspace);
  xi = vector<double>(0, 2*bw-1);

  // Run trafo
  printf("Go to Xi.\n"); 
  Naive_SynthesizeX(coeffs, bw, 0, xi, plm);
  

  theta = vector<double>(0, 2*bw-1);
  ArcCosEvalPts(2*bw, theta);

  //for(l=0; l<2*bw; l++) printf("%g %g\n", theta[l], xi[l]);
  
  // Prepare for trafo back:
  printf("Prepare to go back.\n");
  weights = vector<double>(0,4*bw-1);
  makeweights(bw, weights);
  // Trafo back!
  printf("Transform back\n");
  newcoeffs = vector<double>(0,bw-1);
  Naive_AnalysisX(xi, bw, 0, weights, newcoeffs, plm, workspace);
  // Apply missing factors:
  printf("Apply factors.\n");
  ApplyClFactors(newcoeffs, bw, 8000, 6);
  
  //for(l=0; l<bw; l++) printf("%g\n", newcoeffs[l]);
  
  free_vector(weights, 0, 4*bw-1);
  free_vector(newcoeffs, 0, bw-1);
  free_vector(coeffs, 0, maxL-1);
  free_vector(xi, 0, 2*bw-1);
  free_vector(theta, 0, 2*bw-1);
  free_vector(ll,  0, ClLength-1);
  free_vector(Cov, 0, ClLength-1);
  free_vector(workspace, 0, 16*bw-1);
  free_vector(plm, 0, 2*bw*bw-1);
}

