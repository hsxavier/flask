#ifndef FITSFUNCTIONS_H
#define FITSFUNCTIONS_H 1

#include "ParameterList.hpp"

int WriteCatalog2Fits(std::string filename, double **table, long Nentries, const ParameterList & config);

#endif
