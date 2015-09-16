#ifndef CLPROCESSING     // include guard.
#define CLPROCESSING 1

#include <gsl/gsl_matrix.h>     // gsl_matrix.
#include "FieldsDatabase.hpp"
#include "ParameterList.hpp"

int ClProcess(gsl_matrix ***CovBylAddr, int *Nls, const FZdatabase & fieldslist, const ParameterList & config);

std::string PrintOut(std::string prefix, int i, int j, const FZdatabase & filelist, double *x, double *y, int length);
int GetGaussCorr(double *gXi, double *lnXi, int XiLength, double mean1, double shift1, double mean2, double shift2);
void GetLNCorr(double *lnXi, double *gXi, int XiLength, double mean1, double shift1, double mean2, double shift2);

#endif
