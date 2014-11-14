#include <gsl/gsl_sf_legendre.h>
#include <math.h>
#include "../src/Utilities.hpp"
#include "../src/interpol.h"           // Interpolation. 
#include <stdio.h>

void LegendreC2F(double *C, int maxl, double *F, double *costheta, int ntheta) {
  const double overfourpi = 0.0795774715459;
  double sum=0;
  int l, itheta;

  for (itheta=0; itheta<ntheta; itheta++) {
    F[itheta] = 0;
    for(l=0; l<maxl; l++) F[itheta] += ((double)(2*l+1))*C[l]*gsl_sf_legendre_Pl(l, costheta[itheta]);
    F[itheta] = F[itheta]*overfourpi;
  }
}

double acosCheb(int k, int max){
  const double pi = 3.14159265359;
  return ((double)(2*k-1))/((double)(2*max))*pi;
}

double Cheb(int k, int max) {
  return cos(acosCheb(k, max));
}

int main () {
  double *ll, *Cov, *wrapper[2], *coeffs, *xi, *theta, *costheta; 
  int bw=1000, ClLength=200, maxL=bw, l;
  
  // Load C_l:
  ll  = vector<double>(0, ClLength-1);
  Cov = vector<double>(0, ClLength-1);
  wrapper[0]=ll;
  wrapper[1]=Cov;
  ImportVecs(wrapper, ClLength, 2, "../data/ps_test_grf.dat");
  coeffs = vector<double>(0, maxL-1);
  for(l=0; l<maxL; l++) coeffs[l] = Interpol(ll, ClLength, Cov, (double)l);
  coeffs[0]=0; // Set C_0 = 0
  free_vector(ll,  0, ClLength-1);
  free_vector(Cov, 0, ClLength-1);

  // Set theta:
  theta    = vector<double>(0, 2*bw-1);
  costheta = vector<double>(0, 2*bw-1);
  for (l=0; l<2*bw; l++) {
    theta[l]    = acosCheb(l+1, 2*bw);
    costheta[l] = cos(theta[l]);
  }
  
  // Compute Xi:
  xi = vector<double>(0, 2*bw-1);
  LegendreC2F(coeffs, maxL, xi, costheta, 2*bw);
  
  for(l=0; l<=2*bw; l++) printf("%g %g\n", theta[l], xi[l]);
    
  return 0;
}
