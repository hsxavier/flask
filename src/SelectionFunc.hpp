#ifndef SELECTIONFUNC_H    // include guard.
#define SELECTIONFUNC_H 1

#include "ParameterList.hpp"
#include <healpix_map.h>
#include <gsl/gsl_randist.h>    // Random numbers.

// SelectionFunction class interface:
class SelectionFunction {
private:
  Healpix_Map<double> *AngularSel;
  Healpix_Map<double>  StarMask;
  double **zSel, **zEntries, **fieldZrange, zSearchTol, Scale;
  long *NzEntries;
  int Separable, Nfields, *ftype, Npixels, N1, N2, *zSelIndex, NgalTypes, SelectionType;
  int IndexGalTypes();
public:
  SelectionFunction();
  void load(const ParameterList & config, int *ftype0, double **fzrange, int N10, int N20);
  int Nside();
  int Scheme();
  double operator()(int fz, int pix);
  double operator()(int fz);
  double RandRedshift(gsl_rng *r, int fz, int pix);
  ~SelectionFunction();
};


// Other functions, not members but related:
void SelectionMemTest1(const ParameterList & config, int *ftype0, double **fzrange, int N10, int N20);
#endif
