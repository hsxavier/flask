#ifndef FITSFUNCTIONS_H
#define FITSFUNCTIONS_H 1

#include "definitions.hpp"
#include "ParameterList.hpp"
#include <arr.h>
 
int WriteCatalog2Fits(std::string filename, CAT_PRECISION **table, long Nentries, const ParameterList & config);
int ReadHealpixWeights(int col, int nside, const ParameterList & config, double *weights);

#endif
