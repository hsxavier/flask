#ifndef DEFINITIONS
#define DEFINITIONS 1

// Precision parameters

#define ALM_PRECISION double // Precision for Healpix alm's, float or double. 
#define MAP_PRECISION double // Precision for Healpix maps, float or double.
#define SEL_PRECISION double // Precision for Healpix maps used as mask and selection function.
#define CAT_PRECISION double // Presition for the catalog (matrix), float or double.

// Global definitions

namespace definitions {
  const int fgalaxies=1, fshear=2;                      // Field type identification.
  enum simtype {gaussian, lognormal, homogeneous};      // Type of simulation.
}


#endif
