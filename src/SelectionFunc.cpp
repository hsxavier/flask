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
  NoMap=-2;     // Interval constant.
}

// Load selection functions:
void SelectionFunction::load(const ParameterList & config, int *ftype0, double **fzrange, const FZdatabase & fieldlist) {
  using namespace definitions;
  std::string tempstr, filename, starfile;
  char message[100];
  int i, Nside=-1, scheme=-1, f, z, prevf;
  double *wrapper[2];
  long Ncolumns;

  // Overall properties of the selection functions and fields:
  //N1 = N10;     N2 = N20;
  Nfields       = fieldlist.Nfields();
  zSearchTol    = config.readd("ZSEARCH_TOL");
  SelectionType = config.readi("SELEC_TYPE");
  Separable     = config.readi("SELEC_SEPARABLE");
  tempstr       = config.reads("SELEC_PREFIX");
  starfile      = config.reads("STARMASK");
  Scale         = config.readd("SELEC_SCALE");
  fieldZrange   = matrix<double>(0, Nfields-1, 0, 1);
  ftype         = vector<int>   (0, Nfields-1);
  for (i=0; i<Nfields; i++) {
    fieldZrange[i][0] = fzrange[i][0];
    fieldZrange[i][1] = fzrange[i][1];
    ftype[i] = ftype0[i];
  }
  // No selection options:
  if (tempstr =="0") UseAngularMask=0; else UseAngularMask=1;
  if (starfile=="0") UseStarMask   =0; else UseStarMask   =1;

  // Barriers against unimplemented selection function kinds:
  if (SelectionType==1 || SelectionType==3) error("SelectionFunction.load: SELEC_TYPE fraction of gals not implemented yet.");
  if (SelectionType >3 || SelectionType <0) error("SelectionFunction.load: unknown selection function type.");
  if (Separable     >1 || Separable     <0) error("SelectionFunction.load: unkown SELEC_SEPARABLE option.");
  if (SelectionType==2 && Separable    ==0) error("SelectionFunction.load: bookkeeping for non-separable sel. func. not implemented yet.");

  
  // Selection function f_i(z,theta,phi) not separable: need one map per redshift z per galaxy type i:  
  if(Separable==0) {
    if (UseAngularMask==1) {
      // Read selection functions from FITS files:
      AngularSel = vector<Healpix_Map<SEL_PRECISION> >(0,Nfields-1);
      for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	  // Load FITS file:
	  fieldlist.Index2Name(i, &f, &z);
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
  }
  
  
  // Selection function f_i(z,theta,phi) = f_i(z) * m(theta,phi):
  else if(Separable==1) {
    
    if(UseAngularMask==1) {
      // Load a fixed angular selection function:
      AngularSel = vector<Healpix_Map<SEL_PRECISION> >(0,0);
      read_Healpix_map_from_fits(tempstr, AngularSel[0]);
      Nside  = AngularSel[0].Nside();
      scheme = AngularSel[0].Scheme(); 
    }

    // Prepare for radial selection functions:
    zSelIndex = vector<int>    (0, Nfields-1);
    NgalTypes = IndexGalTypes(fieldlist);
    zSel      = vector<double*>(0, NgalTypes-1); 
    zEntries  = vector<double*>(0, NgalTypes-1);
    NzEntries = vector<long>   (0, NgalTypes-1);
    tempstr   = config.reads("SELEC_Z_PREFIX");
    if (tempstr=="0") error ("SelectionFunction.load: SELEC_Z_PREFIX set to 0.");

    // LOOP over fields, look for type of tracer (galaxy):
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	fieldlist.Index2Name(i, &f, &z);
	// Found new type of galaxy: load radial selection function:
	sprintf(message, "%sf%d.dat", tempstr.c_str(), f);
	filename.assign(message);
	LoadVecs(wrapper, filename, &(NzEntries[zSelIndex[i]]), &Ncolumns,0); // This allocates memory for zEntries[i] and zSel[i].
	zEntries[zSelIndex[i]] = wrapper[0]; zSel[zSelIndex[i]] = wrapper[1];
	if (Ncolumns!=2) error ("SelectionFunction.load: Expected two columns in file "+filename);
      }
  }  
  // Unknown type of selection function:
  else error("SelectionFunction.load: unkown SELEC_SEPARABLE option.");
  

  // If SELEC_TYPE==fraction of galaxies, multiply it by the mean projected density:
  // Not implemented yet.
  //if (SelectionType==1 || SelectionType==3) {} 


  // Load mask over stars:
  if (UseStarMask==1) {
    read_Healpix_map_from_fits(starfile, StarMask);   
    if (UseAngularMask==1) {
      if (StarMask.Nside() != Nside) 
	error("SelectionFunction.load: STARMASK number of pixels is different from the selection function.");
      if (StarMask.Scheme() != scheme) 
	error("SelectionFunction.load: STARMASK and selection function have different pixel orderings.");
    }
  }   

  // Final settings:
  Npixels = 12*Nside*Nside;
}


// Count types of galaxies and index their radial selection functions:
// This function assumes zSelIndex is already allocated.
int SelectionFunction::IndexGalTypes(const FZdatabase & fieldlist) {
  using namespace definitions;
  int i, f, z, NgalTypes=0, *IsSet;

  // For bookkeeping:
  IsSet = vector<int>(0, Nfields-1);
  for (i=0; i<Nfields; i++) IsSet[i]=0;

  // Loop over all fields (small f)
  for (f=0; f<fieldlist.Nfs(); f++) {
    // If is a galaxy, count it and set an index for its selection function:
    if (ftype[fieldlist.fFixedIndex(f,0)]==fgalaxies) {
      NgalTypes++;
      for (z=0; z<fieldlist.Nz4f(f); z++) {
	i = fieldlist.fFixedIndex(f, z);
	zSelIndex[i] = NgalTypes-1;
	IsSet[i]++;
      }
    }
    // If not a galaxy, the index is -1:
    else for (z=0; z<fieldlist.Nz4f(f); z++) {
	i = fieldlist.fFixedIndex(f, z);
	zSelIndex[i] = -1;
	IsSet[i]++;
      }
  }
  
  // Check if every Field has been assigned one and only one zSel (including none=-1):
  for (i=0; i<Nfields; i++) 
    if (IsSet[i]!=1) error("SelectionFunction.IndexGalTypes: radial selection function assignment is messed up.");
  
  free_vector(IsSet, 0, Nfields-1);
  return NgalTypes;
}


// Returns Nside (Healpix) of the angular part of the selection function.
int SelectionFunction::Nside() {
  if  (UseAngularMask==1) return AngularSel[0].Nside();
  else if(UseStarMask==1) return StarMask.Nside();
  else return NoMap;
}


// Returns Scheme (Healpix) of the angular part of the selection function.
int SelectionFunction::Scheme() {
  if  (UseAngularMask==1) return (int)(AngularSel[0].Scheme());
  else if(UseStarMask==1) return (int)(StarMask.Scheme());
  else return NoMap;
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
    if (UseAngularMask==1) {
      free_vector(AngularSel, 0, Nfields-1);
    }
  }
  // Deallocate separable selection functions:
  else if(Separable==1) {
    // Free angular part:
    if (UseAngularMask==1) free_vector(AngularSel, 0, 0);
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


int SelectionFunction::MaskBit(int pix) {
  int bit=0;
  //         ISSUE            BITS
  // Removed by selection:     1
  // Removed by star mask:     2
  // 0 < Selection < 1:        4
  // 0 < Star mask < 1:        8
  
  if (UseStarMask==1) {
    if        (StarMask[pix]==0.0)      bit+=2;
    else if   (StarMask[pix]< 1.0)      bit+=8;
  }
  
  if (UseAngularMask) {
    if (Separable==1) {
      if      (AngularSel[0][pix]==0.0) bit+=1;
      else if (AngularSel[0][pix]< 1.0) bit+=4;
    } 
    else warning("SelectionFunction.MaskBit: not implemented for non-separable selection functions.");
  }

  return bit;
}


// Returns the selection function for the field (or redshift) fz and the angular position pix:
double SelectionFunction::operator()(int fz, int pix) {
  using namespace definitions;
  double z0, zSelInterp, StarValue, AngularValue;

  // Error checks:
  if (ftype[fz]!=fgalaxies) error("SelectionFunction.operator(): this is only set for fields of type galaxies.");
  else if (this->Nside()!=NoMap && (pix >= Npixels || pix < 0)) error("SelectionFunction.operator(): requested pixel is out of range.");
  else if (fz  >= Nfields || fz  < 0) error("SelectionFunction.operator(): unknown requested field.");
  else if (Separable!=0 && Separable!=1) error("SelectionFunction.operator(): this object was not initialized (use load first).");
  
  // If Selection if only used for bookkeeping, do not enforce it here:
  if (SelectionType==2 || SelectionType==3) {
    StarValue    = 1;
    AngularValue = 1;                  // << Note below that non-separables are not affected by this.
  }
  else {
    if (UseStarMask==1)    StarValue    = StarMask[pix];
    else                   StarValue    = 1; 
    if (UseAngularMask==1) AngularValue = AngularSel[0][pix];
    else                   AngularValue = 1;
  }
  
  // Normal execution:
  
  // For non-separable selection functions, return the Healpix map value:
  if (Separable==0) {
    if (UseAngularMask==1) return Scale * StarValue * AngularSel[fz][pix];
    else                   return Scale * StarValue;
  }
  
  // For separable selection functions, multiply radial to angular part:
  else if (Separable==1) {
    z0         = (fieldZrange[fz][0] + fieldZrange[fz][1])/2.0; // Get mean redshift of the field.
    zSelInterp = Interpol(zEntries[zSelIndex[fz]], NzEntries[zSelIndex[fz]], zSel[zSelIndex[fz]], z0);
    if (zSelInterp < 0) error("SelectionFunction.operator(): negative radial selection.");
    return Scale * zSelInterp * StarValue * AngularValue;
  }
}


// Function to test for memory leackage by loading and unloading selection functions:
void SelectionMemTest1(const ParameterList & config, int *ftype0, double **fzrange, const FZdatabase & fieldlist) {
  SelectionFunction test;

  // Load Selection function to allocate memory:
  test.load(config, ftype0, fzrange, fieldlist);
  // Destructor should be called when exiting this function.
}
