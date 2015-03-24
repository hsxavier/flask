#ifndef CORRLNFIELDS_AUX     // include guard.
#define CORRLNFIELDS_AUX 1

#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include "Utilities.hpp"
#include <math.h>
#include <gsl/gsl_randist.h> // for ran_redshift function.
#include "Cosmology.hpp" // for ran_redshift function.
#include <vec3.h>        // For xyz2ang function.
#include <pointing.h>    // For xyz2ang function.
#include <healpix_map.h> // For RandAngInPix function.


// Global definitions
namespace definitions {
  const int fgalaxies=1, fshear=2;                      // Field type identification.
}


// Auxiliary functions for corrlnfields program:
void GenEllip(gsl_rng *r, double sigma, double kappa, double gamma1, double gamma2, double *eps1, double *eps2);
pointing RandAngInPix(gsl_rng *r, const Healpix_Map<double> & map, int pixel);
pointing randang(gsl_rng *r, double thetamin, double thetamax, double phimin, double phimax);
pointing xyz2ang(const vec3 & cartesian);
vec3 VecInRotBasis(const pointing & ang, const vec3 & orig);
double RandRedshift0(gsl_rng *r, double zmin, double zmax);
double ran_redshift(gsl_rng *r, double zmin, double zmax, Cosmology *p);
void CorrGauss(double **gaus1, gsl_matrix *L, double **gaus0);
void GetLNCorr(double *lnXi, double *gXi, int XiLength, double mean1, double shift1, double mean2, double shift2);
int GetGaussCorr(double *gXi, double *lnXi, int XiLength, double mean1, double shift1, double mean2, double shift2);
int GetGaussCov(gsl_matrix *gCovar, gsl_matrix *lnCovar, double *means, double *shifts);
double Gauss2LNvar(double gvar, double mean, double variance, double shift);
int getll(const std::string filename);
std::string getllstr(const std::string filename);
void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2); // Not used by 08-jan-2015.
void CountEntries(std::string filename, long *nr, long *nc); // Not used by 08-jan-2015.
void fz2n (int a1, int a2, int *n, int N1, int N2);
void n2fz (int n, int *a1, int *a2, int N1, int N2);
void fzfz2ij (int a1, int a2, int b1, int b2, int *i, int *j, int N1, int N2); // Not used by 08-jan-2015.
void ij2fzfz (int i, int j, int *a1, int *a2, int *b1, int *b2, int N1, int N2); // Not used by 08-jan-2015.
void test_fzij (int N1, int N2);
std::string PrintOut(std::string prefix, int i, int j, int N1, int N2, double *x, double *y, int length);

#endif
