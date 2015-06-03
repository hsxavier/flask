#ifndef GENERALOUTPUT     // include guard.
#define GENERALOUTPUT 1

#include <healpix_map.h> // For RandAngInPix function.
#include <alm.h>         // For GeneralOutput function.
#include <xcomplex.h>    // For GeneralOutput function.
#include "ParameterList.hpp"
#include "Utilities.hpp"
#include <gsl/gsl_matrix.h> // For GSL matrix output.

/*** Matrix output ***/

// Prints a GSL matrix to a file:
void GeneralOutput(const gsl_matrix *Cov, std::string filename, bool inform = 1);
// Prints all GSL matrices in a vector to files:
void GeneralOutput(gsl_matrix **CovByl, const ParameterList & config, std::string keyword, bool inform = 1);


/*** Cl's output ***/

// Prints all Cl's to a TEXT file:
void GeneralOutput(double **recovCl, int N1, int N2, const ParameterList & config, std::string keyword, bool inform = 1);


/*** Alm's output ***/

// Many alm's in one single table:
void GeneralOutput(Alm<xcomplex <double> > *af, const ParameterList & config, std::string keyword, int N1, int N2, bool inform = 1);
// One alm in one table, named with a prefix and the field ID:
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, int f, int z, bool inform = 1);
// One alm in one table, names with a keyword:
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, bool inform = 1);


/*** Healpix map output ***/

// Prints one single map to FITS file based on a PREFIX and a FIELD ID:
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword, int *fnz, bool inform = 1);
// Prints a list of maps to a single TEXT file:
void GeneralOutput(Healpix_Map<double> *mapf, const ParameterList & config, std::string keyword, int N1, int N2, bool fits=0, bool inform = 1);
// Prints two lists of maps to a single TEXT file:
void GeneralOutput(Healpix_Map<double> *gamma1, Healpix_Map<double> *gamma2, 
		   const ParameterList & config, std::string keyword, int N1, int N2, bool inform = 1);
// One set of (Kappa, gamma1, gamma2) maps to one FITS file, named with a prefix and the field ID:
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, int f, int z, bool inform = 1);
// One set of (Kappa, gamma1, gamma2) maps to one FITS file, named with a keyword:
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, bool inform = 1);
// One single map to one FITS file, named with a keyword:
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword, bool inform = 1);
#endif
