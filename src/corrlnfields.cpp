/* corrlnfields: Written by Henrique S. Xavier on Nov-2014
   e-mail: hsxavier@if.usp.br
 */

#include <iostream>
#include <gsl/gsl_linalg.h>     // Cholesky descomposition.
#include <gsl/gsl_randist.h>    // Random numbers.
#include <cstdlib>              // For function 'popen'.
#include <iomanip>              // For 'setprecision'.
#include <alm.h>
#include <xcomplex.h>
#include <healpix_map.h>
#include <alm_healpix_tools.h>
#include <healpix_map_fitsio.h>
#include <levels_facilities.h>
#include <vec3.h>
#include <gsl/gsl_eigen.h>      // debug
#include "corrlnfields_aux.hpp" // Auxiliary functions made for this program.
#include "GeneralOutput.hpp"    // Various file output functions.
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include "Cosmology.hpp"        // Parameters and formulas.
#include "SelectionFunc.hpp"
#include "RegularizeCov.hpp"
#include "ClProcessing.hpp"

/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl;                     // Basic stuff.
  using namespace definitions;                          // Global definitions.
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  Cosmology cosmo;                                      // Cosmological parameters.
  char message[100];                                    // Handling warnings and errors.
  std::string filename;
  std::ofstream outfile;                                // File for output.
  simtype dist;                                         // For specifying simulation type.
  gsl_matrix **CovByl; 
  int status, i, j, l, m, N1, N2, Nfields, mmax, *ftype, Nls;
  double *means, *shifts, **zrange; 
  long long1, long2;
  gsl_set_error_handler_off();                                              // !!! All GSL return messages MUST be checked !!!
  

  /**********************************************/
  /*** PART 0: Test code and load config file ***/
  /**********************************************/

  // Testing the code:
  cout << "Testing the code... "; cout.flush();
  test_fzij(3,11); test_fzij(13,4);
  cout << "done.\n";
  
  // Loading config file:
  if (argc<=1) { cout << "You must supply a config file." << endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl;
  cosmo.load(&config);
  if (config.reads("EXIT_AT")!="0") config.findpar(config.reads("EXIT_AT")+":"); // Produce warning if last output is unknown.
  // - Lognormal or Gausian realizations:
  if (config.reads("DIST")=="LOGNORMAL") dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN") dist=gaussian;
  else error("corrlnfields: unknown DIST: "+config.reads("DIST"));
  

  /***********************************/
  /*** PART 1: Loads fields info   ***/
  /***********************************/
  bool *fnzSet;
  double **aux;

  // Load means, shifts, type and z range data file:
  cout << "Loading means and shifts from file "+config.reads("FIELDS_INFO")+"... "; cout.flush();
  aux     = LoadTable<double>(config.reads("FIELDS_INFO"), &long1, &long2);
  Nfields = (int)long1;
  fnzSet  = vector<bool>     (0, Nfields-1); for (i=0; i<Nfields; i++) fnzSet[i]=0;
  means   = vector<double>   (0, Nfields-1);
  ftype   = vector<int>      (0, Nfields-1);
  zrange  = matrix<double>   (0, Nfields-1, 0,1);
  if (dist==lognormal) shifts = vector<double>(0, Nfields-1);
  // Check if some input properties are as expected:
  if  (Minimum(aux, 0, Nfields)!=1) error("corrlnfields: first field index in FIELDS_INFO should be 1.");
  if  (Minimum(aux, 1, Nfields)!=1) error("corrlnfields: first redshift index FIELDS_INFO should be 1.");
  N1 = Maximum(aux, 0, Nfields);
  N2 = Maximum(aux, 1, Nfields);
  if  (N1*N2 != Nfields)            error("corrlnfields: mismatch between number and indexes of fields.");
  // Should introduce other checks here, e.g. if there are no gaps or other issues in field IDs.
  
  // Parse information to separate arrays:
  for (j=0; j<Nfields; j++) {
    fz2n((int)aux[j][0], (int)aux[j][1], &i, N1, N2); // Find conventional position of field in arrays.
    if (fnzSet[i]==1) error ("corrlnfields: found more than one mean & shift entry for the same f-z.");
    fnzSet[i] = 1; 
    means[i]  = aux[j][2];  
    if (dist==lognormal) shifts[i] = aux[j][3];
    ftype[i] = (int)aux[j][4];
    zrange[i][0] = aux[j][5]; zrange[i][1] = aux[j][6]; 
  }
  // A few checks on the input:
  for (i=0; i<Nfields; i++) if (fnzSet[i]!=1) error("corrlnfields: the properties of a field were not set.");
  for (i=0; i<Nfields; i++) if (zrange[i][0]>zrange[i][1]) error("corrlnfields: zmin > zmax for a field.");
  if (dist==lognormal) for (i=0; i<Nfields; i++) if(means[i]+shifts[i]<=0) {
	printf(message, "corrlnfields: mean+shift at position %d must be greater than zero.", i); error(message);
      }
  free_vector(fnzSet, 0, Nfields-1);
  free_matrix(aux, 0, Nfields-1, 0, long2-1);
  cout << "done.\n";
  
  cout << "Infered from FIELDS_INFO file:  Nf = " << N1 << "   Nz = " << N2 << endl;
  

  /**************************************************************/
  /*** PART 2: Loads mixing matrices or compute them from Cls ***/
  /**************************************************************/  
  std::string CholeskyInPrefix;
  int lmax, lmin;
  
  CholeskyInPrefix = config.reads("CHOL_IN_PREFIX");
  lmax             = config.readi("LMAX");
  lmin             = config.readi("LMIN");

  // If input triangular mixing matrices unspecified, compute them from input Cls:
  if (CholeskyInPrefix=="0") {
    // Load C(l)s and compute auxiliary Gaussian cov. matrices:
    status = ClProcess(&CovByl, means, shifts, N1, N2, &Nls, config);
    if (status==1) { // Exit if fast output was inside ClProcess.
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0; 
    }
    cout << "Maximum l in input C(l)s: "<<Nls-1<<endl;
    cout << "Will use "<<lmin<<" <= l <= "<<lmax<<endl;
    
    // Cholesky decomposition:
    cout << "Performing Cholesky decompositions of cov. matrices...    "; cout.flush();
    j=0; // Will count number of Cholesky failures.
    for (l=lmin; l<=lmax; l++) {
      //cout << "** Working with cov. matrix for l="<<l<<":\n";
      status = gsl_linalg_cholesky_decomp(CovByl[l]);
      if (status==GSL_EDOM) { 
	sprintf(message,"Cholesky decomposition failed: cov. matrix for l=%d is not positive-definite.", l); 
	warning(message); j++; 
      }
      //cout << "done.\n";
    }
    cout << "done.\n";
    // Exit if any Cholesky failed:
    if (j>0) {sprintf(message,"Cholesky decomposition failed %d times.",j); error(message);}
    // Output mixing matrices if requested:
    GeneralOutput(CovByl, config, "CHOLESKY_PREFIX", 0);
    if (config.reads("CHOLESKY_PREFIX")!="0") 
      cout << ">> Mixing matrices written to prefix "+config.reads("CHOLESKY_PREFIX")<<endl;
  }

  // If input triangular matrices are specified, allocate memory for them:
  else {
    cout << "Allocating memory for mixing matrices (CHOL_IN_PREFIX)... "; cout.flush();
    CovByl = GSLMatrixArray(lmax+1, Nfields, Nfields); // Allocation should have offset to avoid unnecessary low ells.
    cout << "done.\n";                                 // If we are loading the matrices ell by ell, an array is not necessary! 
    cout << "Loading mixing matrices...                                "; cout.flush();
    for (l=lmin; l<=lmax; l++) {
      filename = CholeskyInPrefix+"l"+ZeroPad(l,lmax)+".dat";
      LoadGSLMatrix(filename, CovByl[l]);
    }
    status=0;
    cout << "done.\n";    
  }

  // Exit if dealing with mixing matrices was the last task:
  if (config.reads("EXIT_AT")=="CHOLESKY_PREFIX") {
    cout << "\nTotal number of warnings: " << warning("count") << endl;
    cout<<endl;
    return 0;
  }


  /*************************************************/
  /*** PART 4: Auxiliary Gaussian alm generation ***/
  /*************************************************/
  const double OneOverSqr2=0.7071067811865475;
  bool almout;
  double **gaus0, **gaus1;
  gsl_rng *rnd;
  Alm<xcomplex <double> > *aflm;
  int jmax, jmin;
    
  // Set random number generator:
  rnd = gsl_rng_alloc(gsl_rng_mt19937);
  if (rnd==NULL) error("corrlnfields: gsl_rng_alloc failed!");
  gsl_rng_set(rnd, config.readi("RNDSEED"));    // set random seed

  // Allocate memory for gaussian alm's:
  cout << "Allocating memory for auxiliary gaussian alm's...         "; cout.flush();
  aflm = vector<Alm<xcomplex <double> > >(0,Nfields-1); // Allocate Healpix Alm objects and set their size and initial value.
  for (i=0; i<Nfields; i++) {
    aflm[i].Set(lmax,lmax);
    for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) aflm[i](l,m).Set(0,0);
  }
  cout << "done.\n";

  // LOOP over l's and m's:
  cout << "Generating auxiliary gaussian alm's...                    "; cout.flush(); 
  jmin = (lmin*(lmin+1))/2;
  jmax = (lmax*(lmax+3))/2;
#pragma omp parallel for ordered schedule(dynamic) private(l, m, i, gaus0, gaus1)
  for(j=jmin; j<=jmax; j++) {
    //for (l=lmin; l<=lmax; l++) for (m=0; m<=l; m++) {

    // Allocate temporary memory for random variables:
    gaus0 = matrix<double>(0,Nfields-1, 0,1); // Complex random variables, [0] is real, [1] is imaginary part.
    gaus1 = matrix<double>(0,Nfields-1, 0,1); 

    l = (int)((sqrt(8.0*j+1.0)-1.0)/2.0);
    m = j-(l*(l+1))/2;
    
    // Generate independent 1sigma complex random variables:
#pragma omp ordered
    {
      if (m==0) for (i=0; i<Nfields; i++) {
	  gaus0[i][0] = gsl_ran_gaussian(rnd, 1.0);
	  gaus0[i][1] = 0.0; 
	}                                                      // m=0 are real, so real part gets all the variance.
      else      for (i=0; i<Nfields; i++) {
	  gaus0[i][0] = gsl_ran_gaussian(rnd, OneOverSqr2);
	  gaus0[i][1] = gsl_ran_gaussian(rnd, OneOverSqr2);
	}
    }
    // Generate correlated complex gaussian variables according to CovMatrix:
    CorrGauss(gaus1, CovByl[l], gaus0);
    
    // Save alm to tensor:
    for (i=0; i<Nfields; i++) aflm[i](l,m).Set(gaus1[i][0], gaus1[i][1]);   
    
    // Free temporary memory:
    free_matrix(gaus0,0,Nfields-1,0,1);
    free_matrix(gaus1,0,Nfields-1,0,1);
    
  } // End of LOOP over l's and m's.
  cout << "done.\n";
  free_GSLMatrixArray(CovByl, Nls);
  
  // If requested, write alm's to file:
  GeneralOutput(aflm, config, "AUXALM_OUT", N1, N2);
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="AUXALM_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  /******************************/
  /*** Part 5: Map generation ***/
  /******************************/
  int nside, npixels;
  Healpix_Map<double> *mapf;
  double expmu, gmean, gvar;
  pointing coord;
  char *arg[5];
  char opt1[]="-bar", val1[]="1";

  // Allocate memory for pixel maps:
  cout << "Allocating memory for pixel maps...                       "; cout.flush();
  nside   = config.readi("NSIDE");
  npixels = 12*nside*nside;
  mapf=vector<Healpix_Map<double> >(0,Nfields-1);
  for(i=0; i<Nfields; i++) mapf[i].SetNside(nside, RING); 		
  cout << "done.\n";
  // Generate maps from alm's for each field:
  cout << "Generating maps from alm's...                             "; cout.flush();
  for(i=0; i<Nfields; i++) alm2map(aflm[i],mapf[i]);
  cout << "done.\n";
  // Write auxiliary map to file as a table if requested:
  GeneralOutput(mapf, config, "AUXMAP_OUT", N1, N2);
  
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="AUXMAP_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  // If LOGNORMAL, exponentiate pixels:
  if (dist==lognormal) {
    cout << "LOGNORMAL realizations: exponentiating pixels...          "; cout.flush();
    for (i=0; i<Nfields; i++) {
      gmean = 0; gvar = 0;
      for (j=0; j<npixels; j++) gmean += mapf[i][j];
      gmean = gmean/npixels;
      for (j=0; j<npixels; j++) gvar  += pow(mapf[i][j]-gmean, 2);
      gvar = gvar/(npixels-1);
      expmu=(means[i]+shifts[i])/exp(gvar/2);
      for(j=0; j<npixels; j++) mapf[i][j] = expmu*exp(mapf[i][j])-shifts[i];
    }
    cout << "done.\n";
  }
  // If GAUSSIAN, only add mean:
  else {
    cout << "GAUSSIAN realizations: adding the field mean values to the pixels... "; cout.flush();
    for (i=0; i<Nfields; i++) for(j=0; j<npixels; j++) mapf[i][j] = mapf[i][j] + means[i];
    cout << "done.\n";
  }
  // Free memory for means and shifts:
  if (dist==lognormal) free_vector(shifts, 0, Nfields-1);
  free_vector(means, 0, Nfields-1);
  
  // Write final map to file as a table if requested:
  GeneralOutput(mapf, config, "MAP_OUT", N1, N2);
  // Map output to fits and/or tga files:
  GeneralOutput(mapf, config, "MAPFITS_PREFIX", N1, N2, 1);
  
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="MAP_OUT" ||
      config.reads("EXIT_AT")=="MAPFITS_PREFIX") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  // If requested, compute and write recovered alm's to file:
  if (config.reads("RECOVALM_OUT")!="0") {
    // Clear variables to receive new alm's:
    for (i=0; i<Nfields; i++) for(l=0; l<=lmax; l++) for(m=0; m<=l; m++) aflm[i](l,m).Set(0,0);
    // Compute alm's from map:
    arr<double> weight(2*mapf[0].Nside());
    weight.fill(1);
    cout << "Recovering alm's from map... "; cout.flush();
    for(i=0; i<Nfields; i++) map2alm(mapf[i],aflm[i],weight);
    cout << "done.\n";
    // Output to file:
    GeneralOutput(aflm, config, "RECOVALM_OUT", N1, N2);
    weight.dealloc();
  }
  free_vector(aflm, 0, Nfields-1);
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="RECOVALM_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }
    

  /************************************/
  /*** Part 6: Maps to Observables  ***/
  /************************************/


  /*** Galaxy fields ***/
  
  double PixelSolidAngle=12.56637061435917/npixels; // 4pi/npixels.
  SelectionFunction selection;
  int f, z;
  
  // Read in selection functions from FITS files and possibly text files (for radial part):
  cout << "Reading selection functions from files... "; cout.flush();
  selection.load(config, ftype, zrange, N1, N2); 
  if (selection.Nside()!=mapf[0].Nside()) error("corrlnfields: Selection function and maps have different number of pixels.");
  if (selection.Scheme()!=mapf[0].Scheme()) error("corrlnfields: Selection function and maps have different pixel ordering schemes.");
  cout << "done.\n";

  // Poisson Sampling the galaxy fields:
  if (config.readi("POISSON")==1) {
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	n2fz(i, &f, &z, N1, N2);
	cout << "Poisson sampling f"<<f<<"z"<<z<<"... "; cout.flush();
	for(j=0; j<npixels; j++) mapf[i][j] = gsl_ran_poisson(rnd, selection(i,j)*(1.0+mapf[i][j])/**PixelSolidAngle*/);
	cout << "done.\n";
      }
  }
  // Just generate the expected number density, if requested:
  else if (config.readi("POISSON")==0) {
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	n2fz(i, &f, &z, N1, N2);
	cout << "Using expected number density for f"<<f<<"z"<<z<<"... "; cout.flush();
	for(j=0; j<npixels; j++) mapf[i][j] = selection(i,j)*(1.0+mapf[i][j])/**PixelSolidAngle*/;
	cout << "done.\n";
      }
  }
  else error ("corrlnfields: unknown POISSON option.");
  
  //GeneralOutput(mapf, config, "MAPFITS_PREFIX", fnz, Nfields);

  /*** Lensing fields ***/

  double coeff;
  Healpix_Map<double> *gamma1f, *gamma2f;
  gamma1f = vector<Healpix_Map <double> >(0,Nfields-1);
  gamma2f = vector<Healpix_Map <double> >(0,Nfields-1);
  Alm<xcomplex <double> > Eflm(lmax,lmax), Bflm(lmax,lmax); // Temp memory
  arr<double> weight(2*mapf[0].Nside()); weight.fill(1);    // Temp memory
  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) Bflm(l,m).Set(0,0);     // B-modes are zero for weak lensing.
  
  // LOOP over convergence fields:
  for (i=0; i<Nfields; i++) if (ftype[i]==fshear) {
      n2fz(i, &f, &z, N1, N2);
      cout << "** Will compute shear for f"<<f<<"z"<<z<<":\n";
      
      // Preparing memory:
      cout << "   Allocating and cleaning memory...                    "; cout.flush();
      gamma1f[i].SetNside(nside, RING); gamma1f[i].fill(0);
      gamma2f[i].SetNside(nside, RING); gamma2f[i].fill(0);
      for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) Eflm(l,m).Set(0,0); // E-modes will be set below.
      cout << "done.\n";  
 
      // Get convergence alm's from convergence map:
      cout << "   Transforming convergence map to harmonic space...    "; cout.flush();
      map2alm(mapf[i], Eflm, weight); // Get klm.
      cout << "done.\n";
 
      // Calculate shear alm's from convergence alm's:
      cout << "   Computing shear harmonic coefficients from klm...    "; cout.flush();
      for(l=0; l<2; l++) for (m=0; m<=l; m++) Eflm(l,m).Set(0,0);
      for(l=2; l<=lmax; l++) { // Use Wayne Hu (2000) to get Elm from klm.
	coeff = sqrt( ((double)((l+2)*(l-1))) / ((double)(l*(l+1))) );
	for (m=0; m<=l; m++) Eflm(l,m).Set(coeff*Eflm(l,m).re,coeff*Eflm(l,m).im);
      }
      cout << "done.\n";
      n2fz(i, &f, &z, N1, N2);
      GeneralOutput(Eflm, config, "SHEAR_ALM_PREFIX", f, z);

      // Go from shear E-mode alm's to gamma1 and gamma2 maps:
      cout << "   Transforming harmonic coefficients into shear map... "; cout.flush();
      alm2map_spin(Eflm, Bflm, gamma1f[i], gamma2f[i], 2);
      cout << "done.\n";
      // Write kappa, gamma1 and gamma2 to FITS file:
      n2fz(i, &f, &z, N1, N2);
      GeneralOutput(mapf[i], gamma1f[i], gamma2f[i], config, "SHEAR_FITS_PREFIX", f, z);

    } // End of LOOP over convergence fields.
  weight.dealloc();
  Eflm.Set(0,0); 
  Bflm.Set(0,0);
  // Output shear maps to TEXT tables:
  GeneralOutput(gamma1f, gamma2f, config, "SHEAR_MAP_OUT", N1, N2);

  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="SHEAR_ALM_PREFIX"  ||
      config.reads("EXIT_AT")=="SHEAR_FITS_PREFIX" || 
      config.reads("EXIT_AT")=="SHEAR_MAP_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  /**********************************/
  /*** Part 6: Generate catalog   ***/
  /**********************************/

  double **catalog, esig;
  int Ngalaxies, gali, pixelNgal, PartialNgal, **catSet, ncols;
  pointing ang;
  int ziter, fiter;

  esig = config.readd("ELLIP_SIGMA");
  cout << "Counting galaxies... "; cout.flush();
  Ngalaxies=0; // Find out the total number of galaxies in all fields and redshifts:
  for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) for (j=0; j<npixels; j++) Ngalaxies+=(int)mapf[i][j];
  cout << "done.   # of galaxies: "<<Ngalaxies<<endl;
  
  ncols=10;
  catalog = matrix<double>(0,Ngalaxies-1,0,ncols-1);
  catSet  = matrix<int>(0,Ngalaxies-1,0,ncols-1);
  for (i=0; i<Ngalaxies; i++) for (j=0; j<ncols; j++) catSet[i][j]=0;
  
  // LOOP over 3D bins (pixels and redshifts):
  cout << "Generating catalog... "; cout.flush();
  gali = 0; PartialNgal = 0;
  for (j=0; j<npixels; j++) for(ziter=1; ziter<=N2; ziter++) {
      // Count total number of galaxies of all types in bin:
      pixelNgal = 0; 
      for (fiter=1; fiter<=N1; fiter++) { 
	fz2n(fiter, ziter, &i, N1, N2);
	if (ftype[i]==fgalaxies) pixelNgal+=(int)mapf[i][j];
      }
      
      // LOOP over field IDs:
      for (fiter=1; fiter<=N1; fiter++) {
	fz2n(fiter, ziter, &i, N1, N2);
	// Add entry of type GALAXY:
	if (ftype[i]==fgalaxies) 
	  for(m=0; m<(int)mapf[i][j]; m++) {
	    ang = RandAngInPix(rnd, mapf[i], j);                               // Bookkeeping
	    catalog[gali][0] = ang.theta;                                      catSet[gali][0]++;
	    catalog[gali][1] = ang.phi;                                        catSet[gali][1]++;
	    catalog[gali][2] = RandRedshift0(rnd, zrange[i][0], zrange[i][1]); catSet[gali][2]++;
	    catalog[gali][3] = fiter; /* Field ID (galaxy type) */             catSet[gali][3]++;
	    catalog[gali][9] = j;                                              catSet[gali][9]++;
	    gali++;
	  }
	// Add entry of type SHEAR:
	else if (ftype[i]==fshear) for (m=0; m<pixelNgal; m++) {
	    catalog[PartialNgal+m][4] = mapf[i][j];    catSet[PartialNgal+m][4]++;
	    catalog[PartialNgal+m][5] = gamma1f[i][j]; catSet[PartialNgal+m][5]++;
	    catalog[PartialNgal+m][6] = gamma2f[i][j]; catSet[PartialNgal+m][6]++;
	    GenEllip(rnd, esig, mapf[i][j], gamma1f[i][j], gamma2f[i][j], &(catalog[PartialNgal+m][7]), &(catalog[PartialNgal+m][8]));
	    catSet[PartialNgal+m][7]++; catSet[PartialNgal+m][8]++;
	  }
      } // End of LOOP over field IDs.
 
      PartialNgal+=pixelNgal;
    } // End of LOOP over pixels.  

  // Check if every entry was set once and only once:
  if (gali!=Ngalaxies)        error ("corrlnfields: Galaxy counting is weird (gali != Ngalaxies).");
  if (PartialNgal!=Ngalaxies) error ("corrlnfields: Galaxy counting is weird (PartialNgal != Ngalaxies).");
  for (i=0; i<Ngalaxies; i++) for (j=0; j<ncols; j++) {
      if (catSet[i][j]<1)     error("corrlnfields: Catalog has missing information.");
      if (catSet[i][j]>1)     {cout<<"j: "<<j<<endl; error("corrlnfields: Catalog entry being set more than once.");}
    }
  free_matrix(catSet, 0,Ngalaxies-1,0,ncols-1);
  cout << "done.\n";
  
  // Write catalog to file if requested:
  if (config.reads("CATALOG_OUT")!="0") {
    filename = config.reads("CATALOG_OUT");
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
    else {
      outfile << "# theta, phi, z, fID, convergence, gamma1, gamma2, ellip1, ellip2, pixelID\n";
      PrintTable(catalog, Ngalaxies, ncols, &outfile); 
      outfile.close();
      cout << ">> Catalog written to " << filename << endl;
    }  
  }
  free_matrix(catalog, 0, Ngalaxies-1, 0, ncols-1);
  

  // End of the program
  free_vector(ftype,   0, Nfields-1 );
  free_matrix(zrange,  0, Nfields-1, 0,1);
  free_vector(mapf,    0, Nfields-1 );
  free_vector(gamma1f, 0, Nfields-1 );
  free_vector(gamma2f, 0, Nfields-1 );
  gsl_rng_free(rnd);
  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}
