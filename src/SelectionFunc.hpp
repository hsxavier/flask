#ifndef SELECTIONFUNC_H    // include guard.
#define SELECTIONFUNC_H 1

#include "ParameterList.hpp"
#include <healpix_map.h>

// Parameter list interface:
class SelectionFunction {
private:
  Healpix_Map<double> *AngularSel;
  double **zSel, **zEntries, **fieldZrange;
  long *NzEntries;
  int Separable, Nfields, *ftype, Npixels;
public:
  SelectionFunction();
  void load(const ParameterList & config, int **fnz, int *ftype0, double **fzrange, int Nfields0);
  double operator()(int fz, int pix);
  double operator()(int fz);
  ~SelectionFunction();
};

#endif
