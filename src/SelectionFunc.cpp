#include "SelectionFunc.hpp"
#include "Utilities.hpp"
#include "corrlnfields_aux.hpp" // For definitions namespace and n2fz function.
#include <healpix_map_fitsio.h>
#include "interpol.h"
#include "Maximize.h"

//using std::cout; using std::endl;

// Default constructor:
SelectionFunction::SelectionFunction () {
  Separable=-1; // Kind of selection function not set yet.
}

// Load selection functions:
void SelectionFunction::load(const ParameterList & config, int *ftype0, double **fzrange, int N10, int N20) {
  using namespace definitions;
  std::string tempstr, filename;
  char message[100];
  int i, Nside=-1, scheme=-1, f, z, prevf;
  double *wrapper[2];
  long Ncolumns;

  // Overall properties of the selection functions and fields:
  N1 = N10;     N2 = N20;
  Nfields     = N1*N2;
  zSearchTol  = config.readd("ZSEARCH_TOL");
  Separable   = config.readi("SELEC_SEPARABLE");
  tempstr     = config.reads("SELEC_PREFIX");
  if (Separable==0 || Separable==1) {
    fieldZrange = matrix<double>(0, Nfields-1, 0, 1);
    ftype       = vector<int>   (0, Nfields-1);
    for (i=0; i<Nfields; i++) {
      fieldZrange[i][0] = fzrange[i][0];
      fieldZrange[i][1] = fzrange[i][1];
      ftype[i] = ftype0[i];
    }
  }
  
  // Selection function f_i(z,theta,phi) not separable: need one map per redshift z per galaxy type i:  
  if(Separable==0) {

    // Read selection functions from FITS files:
    AngularSel = vector<Healpix_Map<double> >(0,Nfields-1);
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	// Load FITS file:
	n2fz(i, &f, &z, N1, N2);
	sprintf(message, "%sf%dz%d.fits", tempstr.c_str(), f, z);
	filename.assign(message);
	read_Healpix_map_from_fits(filename, AngularSel[i]);
	// Check if all selection functions have same Nside and Scheme:
	if (Nside==-1) Nside=AngularSel[i].Nside();
	else if (AngularSel[i].Nside() != Nside) 
	  error("SelectionFunction.load: selection FITS files have different number of pixels.");
	if (scheme==-1) scheme=AngularSel[i].Scheme();
	else if (AngularSel[i].Scheme() != scheme) 
	  error("SelectionFunction.load: selection FITS files have different pixel orderings.");
      }
  }

  // Selection function f_i(z,theta,phi) = f_i(z) * m(theta,phi):
  else if(Separable==1) {

    // Load a fixed angular selection function:
    AngularSel = vector<Healpix_Map<double> >(0,0);
    read_Healpix_map_from_fits(tempstr, AngularSel[0]);
    Nside = AngularSel[0].Nside();

    // Prepare for radial selection functions:
    zSelIndex = vector<int>    (0, Nfields-1);
    NgalTypes = IndexGalTypes();
    zSel      = vector<double*>(0, NgalTypes-1); 
    zEntries  = vector<double*>(0, NgalTypes-1);
    NzEntries = vector<long>   (0, NgalTypes-1);
    tempstr = config.reads("SELEC_Z_PREFIX");
    if (tempstr=="0") error ("SelectionFunction.load: SELEC_Z_PREFIX set to 0.");

    // LOOP over fields, look for type of tracer (galaxy):
    prevf=-1;
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	n2fz(i, &f, &z, N1, N2);
	if (f < prevf) error("SelectionFunction.load: n2fz should be monotonic.");
	else if (f > prevf) {
	  // Found new type of galaxy: load radial selection function:
	  sprintf(message, "%sf%d.dat", tempstr.c_str(), f);
	  filename.assign(message);
	  LoadVecs(wrapper, filename, &(NzEntries[zSelIndex[i]]), &Ncolumns,0,1); // This allocates memory for zEntries[i] and zSel[i].
	  zEntries[zSelIndex[i]] = wrapper[0]; zSel[zSelIndex[i]] = wrapper[1];
	  if (Ncolumns!=2) error ("SelectionFunction.load: Expected two columns in file "+filename);
	  prevf = f;
	}
      }
  }  
  // Unknown type of selection function:
  else error("SelectionFunction.load: unkown SELEC_SEPARABLE option.");
  
  // If SELEC_TYPE==FRACTION, multiply it by the mean projected density:
  if (config.reads("SELEC_TYPE")=="FRACTION") error ("SelectionFunction.load: SELEC_TYPE FRACTION not implemented yet.");
  else if (config.reads("SELEC_TYPE")!="DENSITY") error ("SelectionFunction.load: unknown SELEC_TYPE option.");

  // Final settings:
  Npixels = 12*Nside*Nside;
}


// Count types of galaxies and index their radial selection functions:
// This function assumes zSelIndex is already allocated.
int SelectionFunction::IndexGalTypes() {
  using namespace definitions;
  int i, f, z, prevf, NgalTypes=0;

  prevf = -1;
  // Loop over all fields
  for (i=0; i<Nfields; i++) {
    // If is a galaxy, check if it's a new one and set an index for its selection function:
    if (ftype[i]==fgalaxies) {
      n2fz(i, &f, &z, N1, N2);
      if (f < prevf) error("SelectionFunction.CountGalTypes: n2fz should be monotonic.");
      else if (f > prevf) { NgalTypes++; prevf = f; }
      zSelIndex[i] = NgalTypes-1; 
    }
    // If not a galaxy, the index is -1:
    else zSelIndex[i] = -1;
  }

  return NgalTypes;
}


// Returns Nside (Healpix) of the angular part of the selection function.
int SelectionFunction::Nside() {
  return AngularSel[0].Nside();
}


// Returns Scheme (Healpix) of the angular part of the selection function.
int SelectionFunction::Scheme() {
  return (int)(AngularSel[0].Scheme());
}


// Returns random redshift for a field-redshift-bin 'fz' and pixel 'pix' according to the selection function:
// 'zSearchTol' is the tolerance for z position when looking for the Selection function maximum.
double SelectionFunction::RandRedshift(gsl_rng *r, int fz, int pix) {
  double ymax, yguess, zguess, zBinSize;

  zBinSize = fieldZrange[fz][1] - fieldZrange[fz][0];

  if (Separable==0) { // We should interpolate between selection function in the same pixel.
    error("SelectionFunction.RandRedshift: not yet implemented for non-separable selection functions.");
  }
  
  else if (Separable==1) {
    
    ymax = MaxInterp(fieldZrange[fz][0], fieldZrange[fz][1], zSearchTol, 
		     zEntries[zSelIndex[fz]], NzEntries[zSelIndex[fz]], zSel[zSelIndex[fz]]);
    do {
    zguess = fieldZrange[fz][0] + gsl_rng_uniform(r)*zBinSize;
    yguess = gsl_rng_uniform(r)*ymax;
    } while (yguess > Interpol(zEntries[zSelIndex[fz]], NzEntries[zSelIndex[fz]], zSel[zSelIndex[fz]], zguess));

  }

  return zguess;
}


// Destructor:
SelectionFunction::~SelectionFunction() {
  using namespace definitions;
  int i;
  
  // Deallocate non-separable selection functions:
  if(Separable==0) {
    free_vector(AngularSel, 0, Nfields-1);
  }
  // Deallocate separable selection functions:
  else if(Separable==1) {
    // Free angular part:
    free_vector(AngularSel, 0, 0);
    // Free radial selection functions:
    for (i=0; i<NgalTypes; i++) {
      free_vector(zEntries[i], 0, NzEntries[i]-1);
      free_vector(zSel[i],     0, NzEntries[i]-1);
    }
    // Free pointers to functions and counters:
    free_vector(zSelIndex, 0, Nfields-1);
    free_vector(zSel,      0, NgalTypes-1);
    free_vector(zEntries,  0, NgalTypes-1);
    free_vector(NzEntries, 0, NgalTypes-1);
  }
  
  // General deallocation:
  if (Separable==0 || Separable==1) {
    free_matrix(fieldZrange, 0, Nfields-1, 0, 1);
    free_vector(ftype, 0, Nfields-1);
  }

  // If Separable not set, no memory was allocated: do nothing.
}


// Returns the selection function for the field (or redshift) fz and the angular position pix:
double SelectionFunction::operator()(int fz, int pix) {
  using namespace definitions;
  double z0, zSelInterp;
  // Error checks:
  if (ftype[fz]!=fgalaxies) error("SelectionFunction.operator(): this is only set for fields of type galaxies.");
  else if (pix >= Npixels || pix < 0) error("SelectionFunction.operator(): requested pixel is out of range.");
  else if (fz  >= Nfields || fz  < 0) error("SelectionFunction.operator(): unknown requested field.");
  else if (Separable!=0 && Separable!=1) error("SelectionFunction.operator(): this object was not initialized (use load first).");
  // Normal execution:
  else {
    // For non-separable selection functions, return the Healpix map value:
    if (Separable==0) return AngularSel[fz][pix];
    // For separable selection functions, multiply radial to angular part:
    else if (Separable==1) {
      // Get mean redshift of the field:
      z0         = (fieldZrange[fz][0] + fieldZrange[fz][1])/2.0;
      zSelInterp = Interpol(zEntries[zSelIndex[fz]], NzEntries[zSelIndex[fz]], zSel[zSelIndex[fz]], z0);
      if (zSelInterp < 0) error("SelectionFunction.operator(): negative radial selection.");
      return zSelInterp * AngularSel[0][pix];
    }
  }
}


// Function to test for memory leackage by loading and unloading selection functions:
void SelectionMemTest1(const ParameterList & config, int *ftype0, double **fzrange, int N10, int N20) {
  SelectionFunction test;

  // Load Selection function to allocate memory:
  test.load(config, ftype0, fzrange, N10, N20);
  // Destructor should be called when exiting this function.
}
