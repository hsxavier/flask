#include "Utilities.hpp"
#include "ParameterList.hpp"
#include "Cosmology.hpp"
#include "FieldsDatabase.hpp"
#include "definitions.hpp"
#include <omp.h> 
#include <iostream>
#include <unistd.h>             // For access function.
#include "Spline.hpp"

int main(int argc, char *argv[]) {
  using std::cout; using std::endl;                     // Basic stuff.
  using namespace definitions;                          // Global definitions.
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  Cosmology cosmo;                                      // Cosmological parameters.
  FZdatabase fieldlist;
  char message[100];                                    // Handling warnings and errors.
  std::string filename, prefix, outprefix, *filelist;
  std::ofstream outfile;                                // File for output.
  int i, j, k, l, m, Nf, Nz, f, z, Nfields, MaxThreads, allowmiss;
  int NinputCls;
  bool *fnzSet;
  double *tempCl, *wrapper[2];
  // Part 1:
  int af, az, bf, bz, Nlinput, **fnz;
  long *Nentries, ncols;
  double ***ll, ***Cov;
  bool **IsSet;
  time_t start;


  /*************************/
  /*** PART 0: Load data ***/
  /*************************/
  
  start=time(NULL);
  MaxThreads = omp_get_max_threads();
  cout << "Max. # of threads:  "<<MaxThreads<<endl;
  
  // Loading config file:
  if (argc<=2) { cout << "You must supply a config file and a prefix for the output." << endl; return 0;}
  config.load(argv[1]);
  outprefix.assign(argv[argc-1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  cout << "    Output to prefix: "<<outprefix<<endl;
  config.lineload(argc-1, argv);
  config.show();
  cout << endl;
  cosmo.load(config);
 
  // Load means, shifts, type and z range data file:
  Announce("Loading fields information from file "+config.reads("FIELDS_INFO")+"... ");
  fieldlist.Load(config.reads("FIELDS_INFO"));
  Announce();
  Nfields = fieldlist.Nfields();
  Nf      = fieldlist.Nfs();
  Nz      = fieldlist.Nzs();
  cout << "Infered from FIELDS_INFO file:  Nf = " << Nf << "   Nz = " << Nz << endl;  
  k=0; for(f=0; f<Nf; f++) if(fieldlist.ftype(fieldlist.fFixedIndex(f, 0))==fgalaxies) k++;
  cout << "Will use only density fields:   Nf = " << k  << endl;


  /********************************************/
  /*** PART 1: Load C(l)s and organize them ***/
  /********************************************/

  // Get list of the necessary C(l) files:
  prefix    = config.reads("CL_PREFIX"); 
  NinputCls = Nfields*Nfields;
  filelist  = vector<std::string>(0,NinputCls-1);
  Nentries  = vector<long>(0,NinputCls-1);
  Nlinput   = 0;
  // LOOP over all C(l)s:
  for (k=0; k<NinputCls; k++) {
    i = k/Nfields; j = k%Nfields;
    // Get only density Cls:               <<< This is different from main flask code.
    if (fieldlist.ftype(i)==fgalaxies && fieldlist.ftype(j)==fgalaxies) {
      fieldlist.Index2Name(i, &af, &az);
      fieldlist.Index2Name(j, &bf, &bz);
      sprintf(message, "%sf%dz%df%dz%d.dat", prefix.c_str(), af, az, bf, bz);
      // Load file if found:
      if(access(message, R_OK) == 0) {
	filelist[k].assign(message);
	CountEntries(filelist[k], &(Nentries[k]), &ncols);                                   // Get number of Nls.
	if (ncols!=2) error("Dens2KappaCls: wrong number of columns in file "+filelist[k]);
	if (Nentries[k]>Nlinput) Nlinput=Nentries[k];                                        // Record maximum number of ls.
      }
      // Ignore file if it is not found:
      else Nentries[k]=0;
    }
    // Ignore if Cls is not only densities:
    else Nentries[k]=0;
  }
  
  // Allocate memory to store C(l)s:
  // First two indexes are CovMatrix indexes and last is for ll.
  // fnz stores the order that the fields are stored in CovMatrix.
  fnz     =     matrix<int>(0, Nfields-1, 0, 1);                     // Records what field is stored in each element of CovMatrix.
  fnzSet  =    vector<bool>(0, Nfields-1);                           // For bookkeeping.
  ll      = tensor3<double>(0, Nfields-1, 0, Nfields-1, 0, Nlinput); // Records the ll for each C(l) file. 
  Cov     = tensor3<double>(0, Nfields-1, 0, Nfields-1, 0, Nlinput); // Records the C(l) for each C(l) file.
  IsSet   =    matrix<bool>(0, Nfields-1, 0, Nfields-1);             // For bookkeeping.
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) IsSet[i][j]=0;
  for(i=0; i<Nfields; i++) fnzSet[i]=0;
  
  // Read C(l)s and store in data-cube:
  for (k=0; k<NinputCls; k++) if (Nentries[k]>0) {
      i = k/Nfields; j = k%Nfields;
      fieldlist.Index2Name(i, &af, &az);
      fieldlist.Index2Name(j, &bf, &bz);      
      cout << filelist[k] << " goes to ["<<i<<", "<<j<<"]" << endl;
      // Record the order of the fields in CovMatrix:  
      if (fnzSet[i]==0) { fnz[i][0] = af; fnz[i][1] = az; fnzSet[i] = 1; }
      else if (fnz[i][0] != af || fnz[i][1] != az) error("Dens2KappaCls: field order in CovMatrix is messed up!"); 
      if (fnzSet[j]==0) { fnz[j][0] = bf; fnz[j][1] = bz; fnzSet[j] = 1; }
      else if (fnz[j][0] != bf || fnz[j][1] != bz) error("Dens2KappaCls: field order in CovMatrix is messed up!");
      // Import data:
      wrapper[0] =  &(ll[i][j][0]);
      wrapper[1] = &(Cov[i][j][0]);
      ImportVecs(wrapper, Nentries[k], 2, filelist[k].c_str());
      IsSet[i][j]=1; 
    }
  free_vector(Nentries, 0, NinputCls-1);
  free_vector(filelist, 0, NinputCls-1);
  
  // Check if every field was assigned a position in the CovMatrix:           <<< This is different from the main flask code.
  Announce("Sanity checking field allocation in matrix... ");
  for (i=0; i<Nfields; i++) 
    if (fieldlist.ftype(i)==fgalaxies && fnzSet[i]==0) 
      error("Dens2KappaCls: some position in CovMatrix is unclaimed.");
  free_vector(fnzSet, 0, Nfields-1);
  free_matrix(fnz,    0, Nfields-1, 0,1);
  Announce();

  // Error checking and setting remaining Cov elements by simmetry:
  allowmiss = config.readi("ALLOW_MISS_CL");
  if (allowmiss!=1 && allowmiss!=0) warning("Dens2KappaCls: unknown ALLOW_MISS_CL option; will NOT allow missing Cls.");
  if (allowmiss==1) cout << "ALLOW_MISS_CL=1: will set totally missing Cl's to zero.\n";
  Announce("Consistency checking for Cls... ");
  if (IsSet[0][0]!=1) error("Dens2KappaCls: [0,0] is not set.");
  //#pragma omp parallel for schedule(dynamic) private(l, i, j)
  for (k=0; k<NinputCls; k++) {
    i=k/Nfields;  j=k%Nfields;
    // Cl(i,j) was set:
    if (IsSet[i][j]==1) {
      // Check if every Cl is defined for the same ells:
      for (l=0; l<Nlinput; l++) 
	if (ll[i][j][l]!=ll[0][0][l]) { 
	  error("Dens2KappaCls: current implementation requires all Cls to have the same l range.");
	}
      // Check if double set Cls are the same:
      if (IsSet[j][i]==1) {
	for (l=0; l<Nlinput; l++) 
	  if (Cov[i][j][l]!=Cov[j][i][l]) error("Dens2KappaCls: Cl(i,j) and Cl(j,i) are not the same.");    
      }
      // Copy Cl(i,j) to Cl(j,i):
      else { 
	for (l=0; l<Nlinput; l++) { 
	  Cov[j][i][l] = Cov[i][j][l]; 
	  ll [j][i][l] =  ll[i][j][l]; 
	  IsSet[j][i]=1; 
	}
      }
    }
    // Neither Cl(i,j) nor Cl(j,i) was set:
    else if (IsSet[j][i]==0) {
      if (fieldlist.ftype(i)==fgalaxies && fieldlist.ftype(j)==fgalaxies) {
	if (allowmiss==1) {
	  for (l=0; l<Nlinput; l++) { 
	    Cov[i][j][l]=0.0; ll[i][j][l]=ll[0][0][l];
	    Cov[j][i][l]=0.0; ll[j][i][l]=ll[0][0][l];
	  }
	  IsSet[i][j]=1; IsSet[j][i]=1; 
	}
	else error("Dens2KappaCls: missing density Cl.");  
      }
    }
  }
  free_matrix(IsSet, 0, Nfields-1, 0, Nfields-1);
  Announce();

  /*
  Spline densCl;
  int z1, z2; 
  double zmin=1000.0, zmax=0.0, factor, *xinterp, dxinterp, zmean;
  l=1000;
  f=0;
  printf("l=%g\n", ll[0][0][l]);
  outfile.open("origCov.dat");
  for (z1=0; z1<fieldlist.Nz4f(f); z1++) {
    fieldlist.fFixedIndex(f, z1, &i);
    zmean = (fieldlist.zmin(i)+fieldlist.zmax(i))/2.0;
    if (zmean<zmin) zmin=zmean;
    if (zmean>zmax) zmax=zmean;
    for (z2=0; z2<fieldlist.Nz4f(f); z2++) {
      fieldlist.fFixedIndex(f, z2, &j);
      outfile << Cov[i][j][l] << " ";
    }
    outfile << std::endl;
  }
  outfile.close();

  densCl.init(fieldlist, Cov, f, l);
  factor     = 4;
  xinterp    = vector<double>(0, factor*Nfields-1);
  dxinterp   = (zmax-zmin)/((double)(factor*Nfields)-1.0); 
  outfile.open("interpCov.dat");
  for(i=0; i<factor*Nfields; i++) {
    xinterp[i] = zmin + i * dxinterp;
    for(j=0; j<factor*Nfields; j++) {
      xinterp[j] = zmin + j * dxinterp;
      outfile << densCl(xinterp[i], xinterp[j]) << " ";
    }
    outfile << std::endl;
  }
  outfile.close();
  
  wrapper[0]=xinterp;
  outfile.open("interpDom.dat");
  PrintVecs(wrapper, factor*Nfields, 1, &outfile);
  outfile.close();

  return 0;
  */

  /***********************************************************/
  /*** PART 2: Compute convergence Cls through Riemann sum ***/
  /***********************************************************/
  double **KappaWeightTable;
  int zs1, zs2, zl1, zl2, is1, is2, il1, il2;

  // Error checking (density fields must have continuous redshift coverage):
  k = fieldlist.CheckZ4Int();
  if (k==0) error("Dens2KappaCls: no density field found for integrating");
  
  // Compute Kernel:
  Announce("Tabulating integration kernel... ");
  KappaWeightTable = matrix<double>(0, Nfields-1, 0, Nfields-1);
  for (i=0; i<Nfields; i++) 
    for (j=0; j<Nfields; j++) 
      //               s  l
      KappaWeightTable[i][j] = AvgKappaWeightByZ(cosmo, fieldlist.zmin(j), fieldlist.zmax(j), fieldlist.zmax(i)) 
	* (fieldlist.zmax(j)-fieldlist.zmin(j));
  Announce();
  
  /*
  outfile.open("testWeight.dat");
  PrintTable(KappaWeightTable, Nfields, Nfields, &outfile);
  outfile.close();
  return 0;
  */

  tempCl = vector<double>(0, Nlinput-1);
  
  // Get kappa-kappa Cls :
  Announce("Will compute kappa-kappa Cls:");
  printf("\n");
  // LOOP over density fields:
  for(f=0; f<Nf; f++)
    if(fieldlist.ftype(fieldlist.fFixedIndex(f,0))==fgalaxies) {
      // LOOPs over sources redshifts:
      for(zs1=0; zs1<fieldlist.Nz4f(f); zs1++) 
	for(zs2=zs1; zs2<fieldlist.Nz4f(f); zs2++) {
	  is1 = fieldlist.fFixedIndex(f, zs1);
	  is2 = fieldlist.fFixedIndex(f, zs2);
	  // LOOP over multipoles:
#pragma omp parallel for schedule(dynamic) private(zl1, zl2, il1, il2)
	  for(l=0; l<Nlinput; l++) {
	    tempCl[l]=0.0;
	    // LOOPs over lenses redshifts (integration):
	    for(zl1=0; zl1<=zs1; zl1++)
	      for(zl2=0; zl2<=zs2; zl2++) {
		il1 = fieldlist.fFixedIndex(f,zl1);
		il2 = fieldlist.fFixedIndex(f,zl2);
		tempCl[l] += KappaWeightTable[is1][il1]*KappaWeightTable[is2][il2]*Cov[il1][il2][l];
	      } // End over lenses redshifts (integration).
	  } // End over multipoles.
	  // Output to file:
	  fieldlist.fFixedName(f, zs1, &af, &az);
	  fieldlist.fFixedName(f, zs2, &bf, &bz);
	  sprintf(message,"%sf%dz%df%dz%d.dat", outprefix.c_str(), Nf+af, az, Nf+bf, bz);
	  filename.assign(message);
	  if (access(message, R_OK)==0) error("Dens2KappaCls: output file "+filename+" already exists.");
	  outfile.open(message);
	  if (!outfile.is_open()) error("Dens2KappaCls: cannot open file "+filename);
	  wrapper[0]=ll[0][0];
	  wrapper[1]=tempCl;
	  PrintVecs(wrapper, Nlinput, 2, &outfile);
	  outfile.close();
	  cout << "Written Cl to "<<filename<<endl;
	} // End over sources redshifts.
    } // End over galaxy fields.
  Announce();


  // Get density-kappa Cls :
  Announce("Will compute density-kappa Cls:");
  printf("\n");
  // LOOP over density fields:
  for(f=0; f<Nf; f++)
    if(fieldlist.ftype(fieldlist.fFixedIndex(f,0))==fgalaxies) {
      // LOOPs over source and density redshifts:
      for(zs1=0; zs1<fieldlist.Nz4f(f); zs1++)       // << z_dens
	for(zs2=zs1; zs2<fieldlist.Nz4f(f); zs2++) { // << z_kappa
	  is1 = fieldlist.fFixedIndex(f,zs1);
	  is2 = fieldlist.fFixedIndex(f,zs2);
	  // LOOP over multipoles:
#pragma omp parallel for schedule(dynamic) private(zl2, il2)
	  for(l=0; l<Nlinput; l++) {
	    tempCl[l]=0.0;
	    // LOOPs over lens redshift (integration):
	    for(zl2=0; zl2<=zs2; zl2++) {
	      il2 = fieldlist.fFixedIndex(f,zl2);
	      tempCl[l] += KappaWeightTable[is2][il2]*Cov[is1][il2][l];
	    } // End over lenses redshifts (integration).
	  } // End over multipoles.
	  // Output to file:
	  fieldlist.fFixedName(f, zs1, &af, &az);   // dens
	  fieldlist.fFixedName(f, zs2, &bf, &bz);   // kappa
	  sprintf(message,"%sf%dz%df%dz%d.dat", outprefix.c_str(), af, az, Nf+bf, bz);
	  filename.assign(message);
	  if (access(message, R_OK)==0) error("Dens2KappaCls: output file "+filename+" already exists.");
	  outfile.open(message);
	  if (!outfile.is_open()) error("Dens2KappaCls: cannot open file "+filename);
	  wrapper[0]=ll[0][0];
	  wrapper[1]=tempCl;
	  PrintVecs(wrapper, Nlinput, 2, &outfile);
	  outfile.close();
	  cout << "Written Cl to "<<filename<<endl;
	} // End over sources redshifts.
    } // End over galaxy fields.
  Announce();
  
 
  // Memory deallocation:
  free_tensor3(Cov,    0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_tensor3(ll,     0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  PrepareEnd(start);  
  return 0;
}

