#ifndef REGULARIZECOV_H    // include guard.
#define REGULARIZECOV_H 1

#include <gsl/gsl_matrix.h>
#include "ParameterList.hpp"

void RegularizeCov(gsl_matrix * A, const ParameterList & config);

#endif 
