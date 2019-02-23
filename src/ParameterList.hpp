/***************************************************************
2016-08-15: The class 'ParameterList' allows you to create lists of 
parameters to be used in your program and load their values from 
a 'parameter list file'. These parameters must be listed in the 
'ParDef' namespace below. You must set two lists:
1- The names of the parameters as listed in the parameter list 
   file.
2- The type of each parameter, which can be: i1 (int), i2-i5 
   (int[2-5]), d1 (double), d2-d5 (double[2-5]), c (char), 
   s (char[60]) or ph (phrase, sequence of words separated by spaces).
The 'npars' variable must be set to the number of parameters in the 
'par_name' and 'par_type' arrays.   
***************************************************************/

#ifndef PARLIST_H    // include guard.
#define PARLIST_H 1

#define STRSIZE 110

#include <string>
#include <iostream>

namespace ParDef {
  using std::string;
  enum datatype {i1, i2, i3, i4, i5, d1, d2, d3, d4, d5, c, s, ph};
  const string typelabel[13] = {"i1", "i2", "i3", "i4", "i5", "d1", "d2", "d3", "d4", "d5", "c", "s", "ph"};

  // SET HERE THE PARAMETERS OF THE PROGRAM:
  const int     npars=70;
  const string  par_name[npars] = {"RNDSEED", "DIST", "LRANGE", "CROP_CL", "CL_PREFIX", "FLIST_OUT", "AUXALM_OUT", 
				   "SUPPRESS_L", "SUP_INDEX", "XIOUT_PREFIX","FIELDS_INFO", 
				   "GXIOUT_PREFIX", "GCLOUT_PREFIX", "CHOLESKY_PREFIX", "NSIDE", 
				   "MAP_OUT", "RECOVALM_OUT", "AUXMAP_OUT", "LRANGE_OUT", "MMAX_OUT", 
				   "MAPFITS_PREFIX", "FITS2TGA", "OMEGA_m", "OMEGA_L", "W_de", "GALDENSITY",
                                   "SELEC_PREFIX", "SELEC_TYPE", "SELEC_SEPARABLE", "SELEC_Z_PREFIX", 
				   "POISSON", "CATALOG_OUT", "SHEAR_ALM_PREFIX", "SHEAR_FITS_PREFIX", 
				   "SHEAR_MAP_OUT", "ELLIP_SIGMA", "EXIT_AT", "COVL_PREFIX", "EXTRAP_DIPOLE",
				   "REGULARIZE_METHOD", "NEW_EVAL", "REGULARIZE_STEP", "REG_COVL_PREFIX", 
				   "REG_CL_PREFIX", "CHOL_IN_PREFIX", "REG_MAXSTEPS", "ZSEARCH_TOL",
				   "STARMASK", "SELEC_SCALE", "ANGULAR_COORD", "MAPWER_OUT", 
				   "MAPWERFITS_PREFIX", "ELLIPFITS_PREFIX", "ELLIP_MAP_OUT", "RECOVAUXCLS_OUT", "RECOVCLS_OUT", 
				   "ADD_FRAC", "DENS2KAPPA", "USE_HEALPIX_WGTS", "DENS2KAPPA_STAT", "ALLOW_MISS_CL", 
				   "BADCORR_FRAC", "SCALE_CLS", "WINFUNC_SIGMA", "SMOOTH_CL_PREFIX", "APPLY_PIXWIN", 
				   "SHEAR_LMAX", "USE_UNSEEN", "MINDIAG_FRAC", "CATALOG_COLS"};
  const int     par_type[npars] = {i1, s, i2, i1, s, s, s, d1, d1, s, s, s, s, s, i1, s, s, s, i2, i1, s, 
				   i1, d1, d1, d1, d1, s, i1, i1, s, i1, s, s, s, s, d1, s, s, i1, i1, d1, 
				   d1, s, s, s, i1, d1, s, d1, i1, s, s, s, s, s, s, d1, i1, i1, s, i1, d1, d1, 
				   d1, s, i1, i1, i1, d1, ph};
  // END OF PARAMETER SETTINGS.
}


const int MAXPARS=70;    // Maximum number of parameters ParameterList can hold.
const int MAXPARNAME=20; // Maximum size of parameter name.
// Type of data a parameter can hold:
union data {
  double dvec[5];  
  double dnum;  
  long ivec[5];
  long inum;  
  char cvec[STRSIZE];
  char cnum;
};
// Properties of a parameter:
struct parameter {
  char name[MAXPARNAME];
  data value;
  int type;
};


// Parameter list interface:
class ParameterList {
private: 
  parameter list[MAXPARS];
  int parloaded;
public:
  ParameterList ();
  ParameterList (const char *filename);
  //ParameterList (const ParameterList & plist);
  int findpar (std::string word) const;
  void load (const char *filename);
  void lineload (int argc, char *argv[]);
  void show (std::ostream * output = &std::cout) const;
  void copy (int par_index, long *value) const;
  void copy (int par_index, double *value) const;
  void copy (int par_index, char *value) const;
  void copy (std::string par_name, long *value) const;
  void copy (std::string par_name, double *value) const;
  void copy (std::string par_name, char *value) const;
  int readi (int par_index, int vec_pos = -1) const;
  int readi (std::string par_name, int vec_pos = -1) const;
  double readd (int par_index, int vec_pos = -1) const;
  double readd (std::string par_name, int vec_pos = -1) const;
  char readc (int par_index, int pos = -1) const;
  char readc (std::string par_name, int pos = -1) const;
  std::string reads (int par_index) const;
  std::string reads (std::string par_name) const;
};

#endif
