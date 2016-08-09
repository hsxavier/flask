#ifndef CORRLNFIELDS_AUX     // include guard.
#define CORRLNFIELDS_AUX 1

#include "definitions.hpp"
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include "Utilities.hpp"
#include <math.h>
#include <gsl/gsl_randist.h> // for ran_redshift function.
#include "Cosmology.hpp" // for ran_redshift function.
#include <vec3.h>        // For xyz2ang function.
#include <pointing.h>    // For xyz2ang function.
#include <healpix_map.h> // For RandAngInPix function.
#include <alm.h>         // For Kappa2ShearEmode function.
#include <xcomplex.h>    // For Kappa2ShearEmode function.
#include "FieldsDatabase.hpp"
#include <ctime>         // For PrepareEnd.

// Auxiliary functions for corrlnfields program:
void TabulateKappaWeight(double **KappaWeightTable, const Cosmology & cosmo, const FZdatabase & fieldlist);
void ChangeCoord(CAT_PRECISION **catalog, int theta_pos, int phi_pos, long Ngalaxies, int coordtype);
bool ComputeShearQ(const ParameterList & config);
double MapMean(const Healpix_Map<MAP_PRECISION> & map);
double MapVariance(const Healpix_Map<MAP_PRECISION> & map, double mean);
double MapSkewness(const Healpix_Map<MAP_PRECISION> & map, double mean, double variance);
using namespace definitions;
void PrintMapsStats(Healpix_Map<MAP_PRECISION> *mapf, const FZdatabase & fieldlist, simtype dist, std::ostream *output = &std::cout);
void RecoverCls(Alm<xcomplex <ALM_PRECISION> > *bflm, const FZdatabase & fieldlist, std::string clsKey, const ParameterList & config);
void RecoverAlmCls(Healpix_Map<MAP_PRECISION> *mapf, const FZdatabase & fieldlist, 
		   std::string almKey, std::string clsKey, const ParameterList & config);
void PrepRingWeights(int col, arr<double> & weight, const ParameterList & config);
double rad2deg(double rad);
double theta2dec(double theta);
double phi2ra(double phi);
int FileFormat(std::string);
int CountWords(const std::string header);
int GetSubstrPos(const std::string field, const std::string header);
void CatalogFill(CAT_PRECISION **catalog, long row, int column, double value, char **catSet);
void CatalogFill(CAT_PRECISION **catalog, long row, int column, double value);
void Kappa2ShearEmode(Alm<xcomplex <ALM_PRECISION> > &Elm, Alm<xcomplex <ALM_PRECISION> > &Klm);
void GenEllip(gsl_rng *r, double sigma, double kappa, double gamma1, double gamma2, double *eps1, double *eps2);
pointing RandAngInPix(gsl_rng *r, const Healpix_Map<MAP_PRECISION> & map, int pixel);
pointing randang(gsl_rng *r, double thetamin, double thetamax, double phimin, double phimax);
pointing xyz2ang(const vec3 & cartesian);
vec3 VecInRotBasis(const pointing & ang, const vec3 & orig);
double RandRedshift0(gsl_rng *r, double zmin, double zmax);
double ran_redshift(gsl_rng *r, double zmin, double zmax, Cosmology *p);
void CorrGauss(double **gaus1, gsl_matrix *L, double **gaus0);
int getll(const std::string filename);
std::string getllstr(const std::string filename);
void fz2n (int a1, int a2, int *n, int N1, int N2);
void n2fz (int n, int *a1, int *a2, int N1, int N2);
void fzfz2ij (int a1, int a2, int b1, int b2, int *i, int *j, int N1, int N2); // Not used by 08-jan-2015.
void ij2fzfz (int i, int j, int *a1, int *a2, int *b1, int *b2, int N1, int N2); // Not used by 08-jan-2015.
void test_fzij (int N1, int N2);


// Template definitions


// Returns the minimum of matrix[offset...offset+length][index], keeping index fixed.
template<typename type>
type Minimum(type **matrix, int index, int length, int offset=0) {
  int i;
  type min;
  min = matrix[offset][index];
  for(i=offset+1; i<offset+length; i++) if (matrix[i][index]<min) min=matrix[i][index];
  return min;
}

// Returns the maximum of matrix[offset...offset+length][index], keeping index fixed.
template<typename type>
type Maximum(type **matrix, int index, int length, int offset=0) {
  int i;
  type max;
  max = matrix[offset][index];
  for(i=offset+1; i<offset+length; i++) if (matrix[i][index]>max) max=matrix[i][index];
  return max;
}



#endif
