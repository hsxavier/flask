#include "SelectionFunc.hpp"
#include "Utilities.hpp"
#include "corrlnfields_aux.hpp" // For definitions namespace.
#include <healpix_map_fitsio.h>

//using std::cout; using std::endl;

// Default constructor:
SelectionFunction::SelectionFunction () {
  Separable=-1; // Kind of selection function not set yet.
}

// Load selection functions:
void SelectionFunction::load(const ParameterList & config, int **fnz, int *ftype0, double **fzrange, int Nfields0) {
  using namespace definitions;
  std::string tempstr, filename;
  char message[100];
  int i, j;
  double *wrapper[2];
  long Nrows, Ncolumns;

  // Overall properties of the selection functions and fields:
  Nfields     = Nfields0;
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
