#ifndef GSL_AUX_H   // include guard
#define GSL_AUX_H 1

#include <iostream>
#include <gsl/gsl_matrix.h>

// Import a table from file to gsl_matrix format:
gsl_matrix *LoadGSLMatrix(std::string filename);

// Print GSL matrix as a table (with rows and columns):
void PrintGSLMatrix(const gsl_matrix *A, std::ostream *output = &std::cout);

#endif
