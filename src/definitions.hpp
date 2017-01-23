#ifndef DEFINITIONS
#define DEFINITIONS 1

// Precision parameters: Healpix alm's, Healpix maps, Healpix selection function map, catalog (matrix).
#define ALM_PRECISION float
#define MAP_PRECISION float
#define SEL_PRECISION float
#define CAT_PRECISION float
#define FIT_PRECISION TFLOAT

// Global definitions
namespace definitions {
  const MAP_PRECISION UNSEEN = -1.6375e+30;                  // Healpix value for masked regions.
  const int fgalaxies=1, flensing=2;                         // Field type identification.
  enum simtype {gaussian, lognormal, homogeneous};           // Type of simulation.
  const int unknown_format=0, ascii_format=1, fits_format=2; // File format for output catalogue.      
}


#endif
