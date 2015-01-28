#ifndef GENERALOUTPUT     // include guard.
#define GENERALOUTPUT 1

#include <healpix_map.h> // For RandAngInPix function.
#include <alm.h>         // For GeneralOutput function.
#include <xcomplex.h>    // For GeneralOutput function.
#include "ParameterList.hpp"
#include "Utilities.hpp"

void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, int *fnz);
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, int *fnz);
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword);
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword);
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword);

#endif
