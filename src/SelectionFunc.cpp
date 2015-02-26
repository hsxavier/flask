#include "SelectionFunc.hpp"
#include "Utilities.hpp"
#include "corrlnfields_aux.hpp" // For definitions namespace.
#include <healpix_map_fitsio.h>
#include "interpol.h"

//using std::cout; using std::endl;

// Default constructor:
SelectionFunction::SelectionFunction () {
  Separable=-1; // Kind of selection function not set yet.
}

// Load selection functions:
void SelectionFunction::load(const ParameterList & config, int **fnz, int *ftype0, double **fzrange, int N10, int N20) {
  using namespace definitions;
  std::string tempstr, filename;
  char message[100];
  int i, j, Nside=-1;
  double *wrapper[2];
  long Nrows, Ncolumns;

  // Overall properties of the selection functions and fields:
  N1 = N10;     N2 = N20;
  Nfields     = N1*N2;
  Separable   = config.readi("SELEC_SEPARABLE");
  tempstr     = config.reads("SELEC_PREFIX");
  //cout << "just read normal stuff - Nfields: "<<Nfields<<" Separable: "<<Separable<<endl;
  if (Separable==0 || Separable==1) {
    fieldZrange = matrix<double>(0, Nfields-1, 0, 1);
    ftype       = vector<int>   (0, Nfields-1);
    for (i=0; i<Nfields; i++) {
      fieldZrange[i][0] = fzrange[i][0];
      fieldZrange[i][1] = fzrange[i][1];
      ftype[i] = ftype0[i];
    }
    //cout << "Carregou types e ranges.\n";
  }

  // Selection function f_i(z,theta,phi) not separable: need one map per redshift z per galaxy type i:  
  if(Separable==0) {
    //cout << "não é separável!\n";
    // Read selection functions from FITS files:
    AngularSel = vector<Healpix_Map<double> >(0,Nfields-1);
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	//cout << "vai carregar fits para galáxias - field: "<<i<<endl;
	sprintf(message, "%sf%dz%d.fits", tempstr.c_str(), fnz[i][0], fnz[i][1]);
	filename.assign(message);
	read_Healpix_map_from_fits(filename, AngularSel[i]);
	if (Nside==-1) Nside=AngularSel[i].Nside();
	else if (AngularSel[i].Nside() != Nside) 
	  error("SelectionFunction.load: selection FITS files have different number of pixels.");
      }
    // If SELEC_TYPE==FRACTION, multiply it by the mean projected density:
    if (config.reads("SELEC_TYPE")=="FRACTION") error ("SelectionFunction.load: SELEC_TYPE FRACTION not implemented yet.");
    else if (config.reads("SELEC_TYPE")!="DENSITY") error ("SelectionFunction.load: unknown SELEC_TYPE option.");
  }

  // Selection function f_i(z,theta,phi) = f_i(z) * m(theta,phi):
  else if(Separable==1) {
    //cout << "é separável!\n";
    // Load a fixed angular selection function:
    AngularSel = vector<Healpix_Map<double> >(0,0);
    read_Healpix_map_from_fits(tempstr, AngularSel[0]);
    Nside = AngularSel[0].Nside();
    // Prepare for radial selection functions:
    zSel      = vector<double*>(0, Nfields-1);
    zEntries  = vector<double*>(0, Nfields-1);
    NzEntries = vector<long>   (0, Nfields-1);
    tempstr = config.reads("SELEC_Z_PREFIX");
    if (tempstr=="0") error ("SelectionFunction.load: SELEC_Z_PREFIX set to 0.");
    // LOOP over radial selection files:
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	//cout << "vai carregar z selection para field: "<<i<<endl;
	sprintf(message, "%sf%d.dat", tempstr.c_str(), fnz[i][0]);
	filename.assign(message);
	// Load radial selection functions:
	LoadVecs(wrapper, filename, &(NzEntries[i]), &Ncolumns,0,1); // This allocates memory for zEntries[i] and zSel[i].
	zEntries[i] = wrapper[0]; zSel[i] = wrapper[1];
	if (Ncolumns!=2) error ("SelectionFunction.load: Expected two columns in file "+filename);
	//for (j=0; j<NzEntries[i]; j++) cout << zEntries[i][j] << " " << zSel[i][j] << endl; 
      }
  }
  
  else error("SelectionFunction.load: unkown SELEC_SEPARABLE option.");

  // Final settings:
  Npixels = 12*Nside*Nside;
}


// Destructor:
SelectionFunction::~SelectionFunction() {
  using namespace definitions;
  int i;

  //cout << "entrou no destructor\n";
  // Deallocate non-separable selection functions:
  if(Separable==0) {
    free_vector(AngularSel, 0, Nfields-1);
    //cout << "liberou angulares fits (separavel 0)\n";
  }
  // Deallocate separable selection functions:
  else if(Separable==1) {
    //cout << "entrou separável 1\n";
    // Free angular part:
    free_vector(AngularSel, 0, 0);
    //cout << "liberou angular\n";
    // Free radial selection functions:
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	//cout << "vai liberar radial do field: "<<i<<" Num de entradas: "<<NzEntries[i]<<endl;
	free_vector(zEntries[i], 0, NzEntries[i]-1);
	//cout << "passou um\n";
	free_vector(zSel[i],     0, NzEntries[i]-1);
	//cout << "liberou radial do field: "<<i<<endl;
      }
    // Free pointers to functions and counters:
    //cout << "vai liberar ponteiros\n";
    free_vector(zSel,      0, Nfields-1);
    free_vector(zEntries,  0, Nfields-1);
    free_vector(NzEntries, 0, Nfields-1);
    //cout << "liberou ponteiros\n";
  }
  
  // General deallocation:
  if (Separable==0 || Separable==1) {
    //cout << "vai liberar type e range\n";
    free_matrix(fieldZrange, 0, Nfields-1, 0, 1);
    free_vector(ftype, 0, Nfields-1);
    //cout << "liberou type e range\n";
  }

  // If Separable not set, no memory was allocated: do nothing.
}


// Returns the selection function for the field (or redshift) fz and the angular position pix:
double SelectionFunction::operator()(int fz, int pix) {
  using namespace definitions;
  double z0;
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
      z0 = (fieldZrange[fz][0] + fieldZrange[fz][1])/2.0; 
      return Interpol(zEntries[fz], NzEntries[fz], zSel[fz], z0) * AngularSel[0][pix];
    }
  }
}
