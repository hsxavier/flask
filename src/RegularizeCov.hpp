#ifndef REGULARIZECOV_H    // include guard.
#define REGULARIZECOV_H 1

#include <gsl/gsl_matrix.h>
#include "ParameterList.hpp"

// Compute the maximum fractional difference between two matrices, A/B-1 (element-wise):
double MaxFracDiff(gsl_matrix *A, gsl_matrix *B);

void AbsSort(gsl_vector *vec, int lini, int linf);
void RegularizeCov(gsl_matrix * A, const ParameterList & config);

#endif 
