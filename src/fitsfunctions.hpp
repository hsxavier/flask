#ifndef FITSFUNCTIONS_H
#define FITSFUNCTIONS_H 1

#include "definitions.hpp"
#include "ParameterList.hpp"

int WriteCatalog2Fits(std::string filename, CAT_PRECISION **table, long Nentries, const ParameterList & config);

#endif
