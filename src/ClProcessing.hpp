#ifndef CLPROCESSING     // include guard.
#define CLPROCESSING 1

#include <gsl/gsl_matrix.h>     // gsl_matrix.

int ClProcess(gsl_matrix **CovByl, double *means, double *shifts, int N1, int N2, int *Nls, const ParameterList & config);

std::string PrintOut(std::string prefix, int i, int j, int N1, int N2, double *x, double *y, int length);
void CountEntries(std::string filename, long *nr, long *nc);
void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2);
int GetGaussCorr(double *gXi, double *lnXi, int XiLength, double mean1, double shift1, double mean2, double shift2);
void GetLNCorr(double *lnXi, double *gXi, int XiLength, double mean1, double shift1, double mean2, double shift2);

#endif
