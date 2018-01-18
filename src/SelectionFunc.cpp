#include "SelectionFunc.hpp"
#include "Utilities.hpp"
#include "flask_aux.hpp" // For definitions namespace and n2fz function.
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
void SelectionFunction::load(const ParameterList & config, const FZdatabase & fieldlist) {
  using namespace definitions;
  std::string angprefix, zprefix, filename, starfile;
  char message[100];
  int i, j, Nside=-1, scheme=-1, f, z, prevf;
  double *wrapper[2];
  long Ncolumns;

  // Overall properties of the selection functions and fields:
  //N1 = N10;     N2 = N20;
  Nfields       = fieldlist.Nfields();
  zSearchTol    = config.readd("ZSEARCH_TOL");
  SelectionType = config.readi("SELEC_TYPE");
  Separable     = config.readi("SELEC_SEPARABLE");
  angprefix     = config.reads("SELEC_PREFIX");
  starfile      = config.reads("STARMASK");
  Scale         = config.readd("SELEC_SCALE");
  fieldZrange   = matrix<double>(0, Nfields-1, 0, 1);
  ftype         = vector<int>   (0, Nfields-1);
  for (i=0; i<Nfields; i++) {
    fieldZrange[i][0] = fieldlist.zmin(i);
    fieldZrange[i][1] = fieldlist.zmax(i);
    ftype[i] = fieldlist.ftype(i);
  }
  // No selection options:
  if (angprefix =="0") UseAngularMask=0; else UseAngularMask=1;
  if (starfile=="0") UseStarMask   =0; else UseStarMask   =1;

  // Barriers against unimplemented selection function kinds:
  if (Scale < 0) warning ("SelectionFunction.load: SELEC_SCALE < 0 will cause problems when Poisson/Gaussian sampling galaxy fields.");
  if (SelectionType==1 || SelectionType==3) error("SelectionFunction.load: SELEC_TYPE fraction of gals not implemented yet.");
  if (SelectionType >3 || SelectionType <0) error("SelectionFunction.load: unknown selection function type.");
  if (Separable     >2 || Separable     <0) error("SelectionFunction.load: unkown SELEC_SEPARABLE option.");
  if (SelectionType==2 && Separable    ==0) error("SelectionFunction.load: bookkeeping for non-separable sel. func. not implemented yet.");

  // Check if a lensing source selection function is required
  // (i.e. if ellipticity maps were demanded and no galaxy fields were specified)
  // (or if noise maps were demanded and no galaxy fields are specified)
  i=CountGalaxyFields(fieldlist);
  if (i==0 && (config.reads("ELLIP_MAP_OUT")!="0" || config.reads("ELLIPFITS_PREFIX")!="0")) yesShearSel=1;
  if (i==0 && (config.reads("MAPWER_OUT")!="0" || config.reads("MAPWERFITS_PREFIX")!="0")) yesShearSel=1;
  else yesShearSel=0;
  
  // Selection function f_i(z,theta,phi) not separable: need one map per redshift z per galaxy type i:
  if(Separable==0) {
    if (UseAngularMask==1) {
      // Read selection functions from FITS files:
      AngularSel = vector<Healpix_Map<SEL_PRECISION> >(0,Nfields-1);
      for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies || yesShearSel) {
	  // Load FITS file:
	  fieldlist.Index2Name(i, &f, &z);
	  sprintf(message, "%sf%dz%d.fits", angprefix.c_str(), f, z);
	  filename.assign(message);
	  read_Healpix_map_from_fits(filename, AngularSel[i]);
	  // Check if all selection functions have same Nside and Scheme:
	  if (Nside==-1) { Nside=AngularSel[i].Nside(); Npixels=12*Nside*Nside; }
	  else if (AngularSel[i].Nside() != Nside) 
	    error("SelectionFunction.load: selection FITS files have different number of pixels.");
	  if (scheme==-1) scheme=AngularSel[i].Scheme();
	  else if (AngularSel[i].Scheme() != scheme) 
	    error("SelectionFunction.load: selection FITS files have different pixel orderings.");
	  // Set all negative values to zero (this allows for UNSEEN to be used as input):
	  for (j=0; j<Npixels; j++) if (AngularSel[i][j]<0) AngularSel[i][j]=0;
	}
    }
  }
  
  
  // Selection function is separable in radial and angular parts, f_i(z,theta,phi) = f_i(z) * m(theta,phi):
  else if(Separable==1) {
    
    if(UseAngularMask==1) {
      // Load a fixed angular selection function:
      AngularSel = vector<Healpix_Map<SEL_PRECISION> >(0,0);
      read_Healpix_map_from_fits(angprefix, AngularSel[0]);
      Nside   = AngularSel[0].Nside();
      scheme  = AngularSel[0].Scheme();
      Npixels = 12*Nside*Nside; 
      // Set all negative values to zero (this allows for UNSEEN to be used as input):
      for (j=0; j<Npixels; j++) if (AngularSel[0][j]<0) AngularSel[0][j]=0;
    }

    // Prepare for radial selection functions:
    tracerIndex = vector<int>    (0, Nfields-1);
    NgalTypes   = IndexGalTypes(fieldlist);         // << This initializes tracerIndex. 
    zSel        = vector<double*>(0, NgalTypes-1); 
    zEntries    = vector<double*>(0, NgalTypes-1);
    NzEntries   = vector<long>   (0, NgalTypes-1);
    intZsel     = vector<double> (0, Nfields-1);
    zprefix     = config.reads("SELEC_Z_PREFIX");
    if (zprefix=="0") error("SelectionFunction.load: SELEC_Z_PREFIX set to 0.");
    // Mark so far unused vectors to ease bug detection:
    for (i=0; i<Nfields; i++) intZsel[i]=0.0;
    for (i=0; i<NgalTypes; i++) { zSel[i]=NULL; zEntries[i]=NULL; NzEntries[i]=0; }

    // LOOP over fields, look for type of tracer (galaxy):
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies || yesShearSel) {
	fieldlist.Index2Name(i, &f, &z);
	if (NzEntries[tracerIndex[i]] == 0) {
	  // Found new type of galaxy: load radial selection function:
	  sprintf(message, "%sf%d.dat", zprefix.c_str(), f);
	  filename.assign(message);
	  LoadVecs(wrapper, filename, &(NzEntries[tracerIndex[i]]), &Ncolumns,0); // This allocates memory for zEntries[i] and zSel[i].
	  zEntries[tracerIndex[i]] = wrapper[0]; zSel[tracerIndex[i]] = wrapper[1];
	  if (Ncolumns!=2) error ("SelectionFunction.load: Expected two columns in file "+filename);
	}
	// Integrate radial selection function inside the bins:
	intZsel[i] = DiscreteIntegral(zEntries[tracerIndex[i]], zSel[tracerIndex[i]], NzEntries[tracerIndex[i]], 
				      fieldZrange[i][0], fieldZrange[i][1]);
      }
  }  
  
  // Selection function is separable in radial and angular parts, but f_i(z,theta,phi) = f_i(z) * m_i(theta,phi):
  // (that is, we have a different angular part for each tracer)
  else if(Separable==2) {

    // Prepare for radial selection functions:
    tracerIndex = vector<int>    (0, Nfields-1);
    NgalTypes   = IndexGalTypes(fieldlist);         // << This initializes tracerIndex. 
    zSel        = vector<double*>(0, NgalTypes-1); 
    zEntries    = vector<double*>(0, NgalTypes-1);
    NzEntries   = vector<long>   (0, NgalTypes-1);
    intZsel     = vector<double> (0, Nfields-1);
    zprefix     = config.reads("SELEC_Z_PREFIX");
    if (zprefix=="0") error("SelectionFunction.load: SELEC_Z_PREFIX set to 0.");
    // Mark so far unused vectors to ease bug detection:
    for (i=0; i<Nfields; i++) intZsel[i]=0.0;
    for (i=0; i<NgalTypes; i++) { zSel[i]=NULL; zEntries[i]=NULL; NzEntries[i]=0; }
    // Prepare for angular selection functions:
    AngularSel = vector<Healpix_Map<SEL_PRECISION> >(0,NgalTypes-1);

    // LOOP over fields, look for type of tracer (galaxy):
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies || yesShearSel) {
	fieldlist.Index2Name(i, &f, &z);
	if (NzEntries[tracerIndex[i]] == 0) {
	  // Found new type of galaxy: 
	  // load radial selection function:
	  sprintf(message, "%sf%d.dat", zprefix.c_str(), f);
	  filename.assign(message);
	  LoadVecs(wrapper, filename, &(NzEntries[tracerIndex[i]]), &Ncolumns,0); // This allocates memory for zEntries[i] and zSel[i].
	  zEntries[tracerIndex[i]] = wrapper[0]; zSel[tracerIndex[i]] = wrapper[1];
	  if (Ncolumns!=2) error ("SelectionFunction.load: Expected two columns in file "+filename);
	  // load angular selection function:
	  if(UseAngularMask==1) {
	    sprintf(message, "%sf%d.dat", angprefix.c_str(), f);
	    filename.assign(message);
	    read_Healpix_map_from_fits(filename, AngularSel[tracerIndex[i]]);
	    // Check if all angular selection functions have same Nside and Scheme:
	    if (Nside==-1) { Nside=AngularSel[i].Nside(); Npixels=12*Nside*Nside; }
	    else if (AngularSel[tracerIndex[i]].Nside() != Nside) 
	      error("SelectionFunction.load: selection FITS files have different number of pixels.");
	    if (scheme==-1) scheme=AngularSel[tracerIndex[i]].Scheme();
	    else if (AngularSel[tracerIndex[i]].Scheme() != scheme) 
	      error("SelectionFunction.load: selection FITS files have different pixel orderings.");
	    // Set all negative values to zero (this allows for UNSEEN to be used as input):
	    for (j=0; j<Npixels; j++) if (AngularSel[tracerIndex[i]][j]<0) AngularSel[tracerIndex[i]][j]=0;
	  }
	}
	// Integrate radial selection function inside the bins:
	intZsel[i] = DiscreteIntegral(zEntries[tracerIndex[i]], zSel[tracerIndex[i]], NzEntries[tracerIndex[i]], 
				      fieldZrange[i][0], fieldZrange[i][1]);
      }
  }  
  // END of new selection function.

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
    else {
      Nside   = StarMask.Nside();
      Npixels = 12*Nside*Nside;
    }
    // Set all negative values to zero (this allows for UNSEEN to be used as input):
    for (j=0; j<Npixels; j++) if (StarMask[j]<0) StarMask[j]=0;
  }   

}


// Count types of galaxies and index their radial selection functions:
// This function assumes tracerIndex is already allocated.
int SelectionFunction::IndexGalTypes(const FZdatabase & fieldlist) {
  using namespace definitions;
  int i, f, z, NgalTypes=0, *IsSet;

  // For bookkeeping:
  IsSet = vector<int>(0, Nfields-1);
  for (i=0; i<Nfields; i++) IsSet[i]=0;

  // Loop over all fields (small f)
  for (f=0; f<fieldlist.Nfs(); f++) {
    // If is a galaxy, count it and set an index for its selection function:
    if (ftype[fieldlist.fFixedIndex(f,0)]==fgalaxies || yesShearSel) {
      NgalTypes++;
      for (z=0; z<fieldlist.Nz4f(f); z++) {
	i = fieldlist.fFixedIndex(f, z);
	tracerIndex[i] = NgalTypes-1;
	IsSet[i]++;
      }
    }
    // If not a galaxy, the index is -1:
    else for (z=0; z<fieldlist.Nz4f(f); z++) {
	i = fieldlist.fFixedIndex(f, z);
	tracerIndex[i] = -1;
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
int SelectionFunction::Nside() const {
  if  (UseAngularMask==1) return AngularSel[0].Nside();
  else if(UseStarMask==1) return StarMask.Nside();
  else return NoMap;
}


// Returns Scheme (Healpix) of the angular part of the selection function.
int SelectionFunction::Scheme() const {
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
		     zEntries[tracerIndex[fz]], NzEntries[tracerIndex[fz]], zSel[tracerIndex[fz]]);
    do {
    zguess = fieldZrange[fz][0] + gsl_rng_uniform(r)*zBinSize;
    yguess = gsl_rng_uniform(r)*ymax;
    } while (yguess > Interpol(zEntries[tracerIndex[fz]], NzEntries[tracerIndex[fz]], zSel[tracerIndex[fz]], zguess));

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
  else if(Separable==1 || Separable==2) {
    // Free angular part:
    if (UseAngularMask==1 && Separable==1) free_vector(AngularSel, 0, 0);
    if (UseAngularMask==1 && Separable==2) free_vector(AngularSel, 0, NgalTypes-1);
    // Free radial selection functions:
    for (i=0; i<NgalTypes; i++) {
      free_vector(zEntries[i], 0, NzEntries[i]-1);
      free_vector(zSel[i],     0, NzEntries[i]-1);
    }
    // Free integral of the radial selection function in the z bin:
    free_vector(intZsel,   0, Nfields-1);
    // Free pointers to functions and counters:
    free_vector(tracerIndex, 0, Nfields-1);
    free_vector(zSel,        0, NgalTypes-1);
    free_vector(zEntries,    0, NgalTypes-1);
    free_vector(NzEntries,   0, NgalTypes-1);
  }
  
  
  // General deallocation:
  if (Separable==0 || Separable==1) {
    free_matrix(fieldZrange, 0, Nfields-1, 0, 1);
    free_vector(ftype, 0, Nfields-1);
  }

  // If Separable not set, no memory was allocated: do nothing.
}


int SelectionFunction::MaskBit(int fz, int pix) const {
  int bit=0;
  //         ISSUE            BITS
  // Removed by selection:     1
  // Removed by star mask:     2
  // 0 < Selection < 1:        4
  // 0 < Star mask < 1:        8
  // Selection > 1:           16
  // Star mask > 1:           32
  
  // Star mask bits:
  if (UseStarMask==1) {
    if        (StarMask[pix]==0.0)      bit+=2;
    else if   (StarMask[pix]< 1.0)      bit+=8;
    else if   (StarMask[pix]> 1.0)      bit+=32;
  }
  
  // Selection function bits:
  if (UseAngularMask) {
    // Separable (but a single angular part):
    if (Separable==1) {
      if      (AngularSel[0][pix]==0.0) bit+=1;
      else if (AngularSel[0][pix]< 1.0) bit+=4;
      else if (AngularSel[0][pix]> 1.0) bit+=16;
    }
    // Separable (and one angular part for each tracer):
    else if (Separable==2) {
      if      (AngularSel[tracerIndex[fz]][pix]==0.0) bit+=1;
      else if (AngularSel[tracerIndex[fz]][pix]< 1.0) bit+=4;
      else if (AngularSel[tracerIndex[fz]][pix]> 1.0) bit+=16;      
    }
    // Non-separable:
    else {
      if      (AngularSel[fz][pix]==0.0) bit+=1;
    }
  }

  return bit;
}


// Returns the selection function for the field (or redshift) fz and the angular position pix:
double SelectionFunction::operator()(int fz, int pix) {
  using namespace definitions;
  double z0, zSelInterp, StarValue, AngularValue;

  // Error checks:
  if (ftype[fz]!=fgalaxies && !yesShearSel) 
    error("SelectionFunction.operator(): this is only set for fields of type galaxies or if ellipticity maps are required.");
  else if (this->Nside()!=NoMap && (pix >= Npixels || pix < 0)) error("SelectionFunction.operator(): requested pixel is out of range.");
  else if (fz  >= Nfields || fz  < 0) error("SelectionFunction.operator(): unknown requested field.");
  else if (Separable!=0 && Separable!=1 && Separable!=2) 
    error("SelectionFunction.operator(): this object was not initialized (use load first).");
  
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
  
  // For separable selection functions (and just one angular part), multiply integral of radial to angular part:
  else if (Separable==1) return Scale * intZsel[fz] * StarValue * AngularValue;
  // For separable selection functions (and a different angular part for each tracer), 
  // multiply integral of radial to the specific angular part:
  else if (Separable==2) return Scale * intZsel[fz] * StarValue * AngularSel[tracerIndex[fz]][pix];
}


// Function to test for memory leackage by loading and unloading selection functions:
void SelectionMemTest1(const ParameterList & config, const FZdatabase & fieldlist) {
  SelectionFunction test;

  // Load Selection function to allocate memory:
  test.load(config, fieldlist);
  // Destructor should be called when exiting this function.
}
