#ifndef GSL_AUX_H   // include guard
#define GSL_AUX_H 1

#include <iostream>
#include <gsl/gsl_matrix.h>

// Allocate an array [0...Nmatrices] of Nrows by Ncols gsl_matrices.
gsl_matrix **GSLMatrixArray(int Nmatrices, int Nrows, int Ncols);

// Free memory allocated with function above.
void free_GSLMatrixArray(gsl_matrix **array, int Nmatrices);

// Import a table from file to gsl_matrix format:
gsl_matrix *LoadGSLMatrix(std::string filename);

// Import a table from file to gsl_matrix format (already allocated):
void LoadGSLMatrix(std::string filename, gsl_matrix *matrix);

// Print GSL matrix as a table (with rows and columns):
void PrintGSLMatrix(const gsl_matrix *A, std::ostream *output = &std::cout);

#endif
