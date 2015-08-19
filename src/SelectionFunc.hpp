#ifndef SELECTIONFUNC_H    // include guard.
#define SELECTIONFUNC_H 1

#include "definitions.hpp"
#include "ParameterList.hpp"
#include <healpix_map.h>
#include <gsl/gsl_randist.h>    // Random numbers.
#include "FieldsDatabase.hpp"

// SelectionFunction class interface:
class SelectionFunction {
private:
  int NoMap;
  Healpix_Map<SEL_PRECISION> *AngularSel;
  Healpix_Map<SEL_PRECISION>  StarMask;
  double **zSel, **zEntries, **fieldZrange, zSearchTol, Scale;
  long *NzEntries;
  int Separable, Nfields, *ftype, Npixels, N1, N2, *zSelIndex, NgalTypes, SelectionType, UseStarMask, UseAngularMask;
  int IndexGalTypes(const FZdatabase & fieldlist);
public:
  SelectionFunction();
  void load(const ParameterList & config, int *ftype0, double **fzrange, const FZdatabase & fieldlist);
  int Nside();
  int Scheme();
  double operator()(int fz, int pix);
  double operator()(int fz);
  double RandRedshift(gsl_rng *r, int fz, int pix);
  int MaskBit(int pix);
  ~SelectionFunction();
};


// Other functions, not members but related:
void SelectionMemTest1(const ParameterList & config, int *ftype0, double **fzrange, int N10, int N20);
#endif
