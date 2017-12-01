/* flask: Written by Henrique S. Xavier on Nov-2014
   e-mail: hsxavier@if.usp.br
   Related paper: Xavier et al. 2016, MNRAS Vol. 459, p. 3693 (arXiv:1602.08503)
 */

#include <iostream>
#include <gsl/gsl_linalg.h>     // Cholesky descomposition.
#include <gsl/gsl_randist.h>    // Random numbers.
#include <iomanip>              // For 'setprecision'.
#include <alm.h>
#include <healpix_map.h>
#include <alm_healpix_tools.h>
#include <omp.h>                // For OpenMP functions, not pragmas.
#include <limits.h>             // For finding out max. value of INT variables.
#include "definitions.hpp"      // Global variables and #defines.
#include "flask_aux.hpp"        // Auxiliary functions made for this program.
#include "GeneralOutput.hpp"    // Various file output functions.
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include "SelectionFunc.hpp"
#include "ClProcessing.hpp"
#include "fitsfunctions.hpp"    // For WriteCatalog2Fits function.
#include "lognormal.hpp"
#include "FieldsDatabase.hpp"
#include <unistd.h> // debugging
#include <ctime> // For timing full code run with StartAll.

#define RAND_OFFSET 10000000  // For generating random numbers in parallel, add multiples of this to seed.


/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl;                     // Basic stuff.
  using namespace definitions;                          // Global definitions.
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  Cosmology cosmo;                                      // Cosmological parameters.
  FZdatabase fieldlist;
  char message[100];                                    // Handling warnings and errors.
  time_t StartAll;                                      // For timing the code run.
  double TotalTime;                                     // For timing the code run.
  bool yesShear;
  std::string filename, ExitAt;
  std::ofstream outfile;                                // File for output.
  simtype dist;                                         // For specifying simulation type.
  gsl_matrix **CovByl; 
  int status, i, j, k, l, m, Nf, Nz, f, z, Nfields, Nls, MaxThreads;
  long long1, long2;
  double *expmu, gvar, gvarl, esig;
  gsl_set_error_handler_off();          // !!! All GSL return messages MUST be checked !!!

  
  /**********************************************/
  /*** PART 0: Test code and load config file ***/
  /**********************************************/

  // Testing the code:
  StartAll = time(NULL);
  Announce("Testing the code... "); 
  // Verify that max. value for INT is not smaller than expected:
  sprintf(message, "%d", INT_MAX); filename.assign(message);
  if (filename.size() < 10) 
    warning("flask: INT_MAX is smaller than expected, may mess parallel random number generator.");
  Announce();
  
  MaxThreads = omp_get_max_threads();
  cout << "FLASK commit:       " << FLASKCOMMIT << endl;
  cout << "Max. # of threads:  " << MaxThreads  << endl;
  if (MaxThreads>210) warning("flask: # of threads too big, may mess parallel random number generator.");

  // Loading config file:
  if (argc<=1) { cout << "You must supply a config file." << endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl;
  cosmo.load(config);
  ExitAt = config.reads("EXIT_AT");
  if (ExitAt!="0") config.findpar(ExitAt+":"); // Produce warning if last output is unknown.
  // - Lognormal or Gausian realizations:
  if (config.reads("DIST")=="LOGNORMAL")        dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN")    dist=gaussian;
  else if (config.reads("DIST")=="HOMOGENEOUS") dist=homogeneous;
  else error("flask: unknown DIST: "+config.reads("DIST"));
 

  /***********************************/
  /*** PART 1: Loads fields info   ***/
  /***********************************/
  
  // Load means, shifts, type and z range data file:
  Announce("Loading fields information from file "+config.reads("FIELDS_INFO")+"... ");
  fieldlist.Load(config.reads("FIELDS_INFO"));
  Announce();
  Nfields = fieldlist.Nfields();
  Nf      = fieldlist.Nfs();
  Nz      = fieldlist.Nzs();
  cout << "Infered from FIELDS_INFO file:  Nf = " << Nf << "   Nz = " << Nz << endl;

 
  
  /**************************************************************/
  /*** PART 2: Loads mixing matrices or compute them from Cls ***/
  /**************************************************************/  
  std::string CholeskyInPrefix;
  int lmax, lmin;
  
  CholeskyInPrefix = config.reads("CHOL_IN_PREFIX");
  lmax             = config.readi("LRANGE", 1);
  lmin             = config.readi("LRANGE", 0);
  if (lmin > lmax) error("LRANGE set in the wrong order.");

  // Skip mixing matrices if generating homogeneous uncorrelated fields: matrices would be zero:
  if (dist!=homogeneous) {
  
    // If input triangular mixing matrices unspecified, compute them from input Cls:
    if (CholeskyInPrefix=="0") {
      // Load C(l)s and compute auxiliary Gaussian cov. matrices:
      status = ClProcess(&CovByl, &Nls, fieldlist, config);
      if (status==1) { // Exit if fast output was inside ClProcess.
	PrepareEnd(StartAll); return 0; 
      }
      cout << "Maximum l in input C(l)s: "<<Nls-1<<endl;
      if (lmax>Nls-1) {
	lmax=Nls-1;
	warning("flask: requested LRANGE upper bound is beyond input data, will use existing data instead.");
      }
      cout << "Will use "<<lmin<<" <= l <= "<<lmax<<endl;

      // Cholesky decomposition:
      Announce("Performing Cholesky decompositions of cov. matrices... ");
      j=0; // Will count number of Cholesky failures.
      for (l=lmin; l<=lmax; l++) {
	//cout << "** Working with cov. matrix for l="<<l<<":\n";
	status = gsl_linalg_cholesky_decomp(CovByl[l]);
	if (status==GSL_EDOM) { 
	  sprintf(message,"Cholesky decomposition failed: cov. matrix for l=%d is not positive-definite.", l); 
	  warning(message); j++; 
	}
      }
      Announce();
    
      // Exit if any Cholesky failed:
      if (j>0) {sprintf(message,"Cholesky decomposition failed %d times.",j); error(message);}
      // Output mixing matrices if requested:
      GeneralOutput(CovByl, config, "CHOLESKY_PREFIX", 0);
      if (config.reads("CHOLESKY_PREFIX")!="0") 
	cout << ">> Mixing matrices written to prefix "+config.reads("CHOLESKY_PREFIX")<<endl;
    }

    // If input triangular matrices are specified, allocate memory for them:
    else {
      Announce("Allocating memory for mixing matrices (CHOL_IN_PREFIX)... ");
      CovByl = GSLMatrixArray(lmax+1, Nfields, Nfields); // Allocation should have offset to avoid unnecessary low ells.
      Announce();                                        // If we are loading the matrices ell by ell, an array is not necessary! 
      Announce("Loading mixing matrices... ");
      for (l=lmin; l<=lmax; l++) {
	filename = CholeskyInPrefix+"l"+ZeroPad(l,lmax)+".dat";
	LoadGSLMatrix(filename, CovByl[l]);
      }
      status=0;
      Announce();    
    }
  } // End of IF not homogeneous.
  else cout << "HOMOGENEOUS realizations: skipped mixing matrix preparation.\n";
 
  // Exit if dealing with mixing matrices was the last task:
  if (ExitAt=="CHOLESKY_PREFIX") {
    PrepareEnd(StartAll); return 0;
  }


  // If necessary, compute the Gaussian field variance from theory:
  if (dist==lognormal) {
    Announce("LOGNORMAL realizations: computing exponential factor... ");
    expmu = vector<double>(0, Nfields-1);
    // LOOP over Fields:
    for (i=0; i<Nfields; i++) {
      gvar = 0;
      // Compute Gaussian field variance from mixing matrices:
      for (l=lmin; l<=lmax; l++) {
	gvarl = 0;
	// Compute Cov(l)_ii from mixing matrices:
	for (j=0; j<=i; j++) gvarl += pow(CovByl[l]->data[i*Nfields+j], 2);
	// Add up multipole contributions:
	gvar += (2.0*((double)l)+1.0)*gvarl;
      }
      gvar /= 4.0*M_PI;
      // Compute the exponential factor:
      expmu[i] = (fieldlist.mean(i)+fieldlist.shift(i))/exp(gvar/2.0);
    } // End of LOOP over Fields.
    Announce();
  }
  


  /*************************************************/
  /*** PART 4: Auxiliary Gaussian alm generation ***/
  /*************************************************/
  const double OneOverSqr2=0.7071067811865475;
  bool almout;
  double ***gaus0, ***gaus1;
  gsl_rng **rnd;
  Alm<xcomplex <ALM_PRECISION> > *aflm, *bflm;
  int jmax, jmin, rndseed0;
    
  // Set random number generators for each thread, plus one for serial stuff [0]:
  // Method is meant to:
  //                    (1) generate aux. alm's fast (in parallel);
  //                    (2) give independent samples for different RNDSEED (parallel seeds may never overlap);
  //                    (3) maintain reproducibility (seeds used in each part of computations must be the same for fixed # of threads).
  Announce("Initializing random number generators... ");
  rndseed0 = config.readi("RNDSEED");
  rnd      = vector<gsl_rng*>(0,MaxThreads+1);
  if (rndseed0 > RAND_OFFSET-1) warning("flask: RNDSEED exceeds RAND_OFFSET-1 in code.");
  for (i=0; i<=MaxThreads; i++) {
    rnd[i] = gsl_rng_alloc(gsl_rng_mt19937);
    if (rnd==NULL) error("flask: gsl_rng_alloc failed!");
    gsl_rng_set(rnd[i], i*RAND_OFFSET+rndseed0);    // set random seed
  }
  Announce();

  // Skip alm generation if creating homogeneous uncorrelated fields: all would be zero:
  if (dist!=homogeneous) {
  
    // Allocate memory for gaussian alm's:
    Announce("Allocating memory for auxiliary gaussian alm's... ");
    gaus0 = tensor3<double>(1,MaxThreads, 0,Nfields-1, 0,1); // Complex random variables, [0] is real, [1] is imaginary part.
    gaus1 = tensor3<double>(1,MaxThreads, 0,Nfields-1, 0,1);  
    aflm  = vector<Alm<xcomplex <ALM_PRECISION> > >(0,Nfields-1);   // Allocate Healpix Alm objects and set their size and initial value.
    for (i=0; i<Nfields; i++) {
      aflm[i].Set(lmax,lmax);
      for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) {
#if USEXCOMPLEX // For compatibility with Healpix versions <=3.20 and >=v3.30.          
	  aflm[i](l,m).Set(0,0);
#else
	  aflm[i](l,m).real(0); 
	  aflm[i](l,m).imag(0);
#endif
	}
    }
    Announce();

    // LOOP over l's and m's together:
    Announce("Generating auxiliary gaussian alm's... ");
    jmin = (lmin*(lmin+1))/2;
    jmax = (lmax*(lmax+3))/2;
#pragma omp parallel for schedule(static) private(l, m, i, k)
    for(j=jmin; j<=jmax; j++) {
    
      // Find out which random generator to use:
      k = omp_get_thread_num()+1;
      // Find out which multipole to compute:    
      l = (int)((sqrt(8.0*j+1.0)-1.0)/2.0);
      m = j-(l*(l+1))/2;
    
      // Generate independent 1sigma complex random variables:
      if (m==0) for (i=0; i<Nfields; i++) {
	  gaus0[k][i][0] = gsl_ran_gaussian(rnd[k], 1.0);
	  gaus0[k][i][1] = 0.0;
	}                                                      // m=0 are real, so real part gets all the variance.
      else      for (i=0; i<Nfields; i++) {
	  gaus0[k][i][0] = gsl_ran_gaussian(rnd[k], OneOverSqr2);
	  gaus0[k][i][1] = gsl_ran_gaussian(rnd[k], OneOverSqr2);
	}
    
      // Generate correlated complex gaussian variables according to CovMatrix:
      CorrGauss(gaus1[k], CovByl[l], gaus0[k]);
  
      // Save alm to tensor:
      for (i=0; i<Nfields; i++) {
#if USEXCOMPLEX // For compatibility with Healpix versions <=3.20 and >=v3.30.          
	aflm[i](l,m).Set(gaus1[k][i][0], gaus1[k][i][1]);
#else
	aflm[i](l,m).real(gaus1[k][i][0]);
	aflm[i](l,m).imag(gaus1[k][i][1]);
#endif
      }   
      
    } // End of LOOP over l's and m's.
    Announce();
    free_GSLMatrixArray(CovByl, Nls);
    free_tensor3(gaus0, 1,MaxThreads, 0,Nfields-1, 0,1);
    free_tensor3(gaus1, 1,MaxThreads, 0,Nfields-1, 0,1);
    // If requested, write alm's to file:
    GeneralOutput(aflm, config, "AUXALM_OUT", fieldlist);
  } // End of IF not homogeneous.

  else cout << "HOMOGENEOUS realizations: skipped alm generation.\n";

  // Exit if this is the last output requested:
  if (ExitAt=="AUXALM_OUT") {
    PrepareEnd(StartAll); return 0;
  }
  // If requested, recover Cls from auxiliary alm's:
  RecoverCls(aflm, fieldlist, "RECOVAUXCLS_OUT", config);
  // Exit if this is the last output requested:
  if (ExitAt=="RECOVAUXCLS_OUT") {
    PrepareEnd(StartAll); return 0;
  }


  /******************************/
  /*** Part 5: Map generation ***/
  /******************************/
  int nside, npixels;
  Healpix_Map<MAP_PRECISION> *mapf;
  pointing coord;
  char *arg[5];
  char opt1[]="-bar", val1[]="1";

  /*** Part 5.1: Generate maps from auxiliar alm's ***/

  // Allocate memory for pixel maps:
  Announce("Allocating memory for pixel maps... "); 
  nside    = config.readi("NSIDE");
  yesShear = ComputeShearQ(config);
  if (nside>sqrt(INT_MAX/12)) warning("flask: NSIDE too large, number of pixels will overflow INT variables");
  npixels = 12*nside*nside;
  mapf=vector<Healpix_Map<MAP_PRECISION> >(0,Nfields-1);
  for(i=0; i<Nfields; i++) mapf[i].SetNside(nside, RING); 		
  Announce();

  // Generate maps from alm's for each field if not creating homogeneous uncorrelated fields:
  if (dist!=homogeneous) {
    Announce("Generating maps from alm's... ");
    for(i=0; i<Nfields; i++) {
      alm2map(aflm[i], mapf[i]);
      // Free alm's after making map if they are not needed (see below at **):
      if (dist==lognormal || (dist==gaussian && yesShear==0)) aflm[i].Set(0,0); 
    }
    Announce();
  }
  // Generate mean maps if creating homogeneous fields:
  else {
    Announce("HOMOGENEOUS realizations: filing maps with mean values... ");
    for(i=0; i<Nfields; i++) mapf[i].fill(fieldlist.mean(i));
    Announce();
  }
  
  // ** If generating lognormal, alm's are not needed anymore (for gaussian the klm's are used to generate shear):
  if (dist==lognormal || (dist==gaussian && yesShear==0)) free_vector(aflm, 0, Nfields-1);
  // Write auxiliary map to file as a table if requested:
  GeneralOutput(mapf, config, "AUXMAP_OUT", fieldlist);
  
  // Exit if this is the last output requested:
  if (ExitAt=="AUXMAP_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  // If LOGNORMAL, exponentiate pixels:
  if (dist==lognormal) {
    Announce("LOGNORMAL realizations: exponentiating pixels... ");
    long2 = ((long)Nfields)*((long)npixels);
#pragma omp parallel for private(i, j)
    for (long1=0; long1<long2; long1++) {
      i = (int)(long1%Nfields);
      j = (int)(long1/Nfields);
      mapf[i][j] = expmu[i]*exp(mapf[i][j])-fieldlist.shift(i);
    }
    free_vector(expmu, 0, Nfields-1);
    Announce();
  }
  

  // If GAUSSIAN, only add mean:
  else if (dist==gaussian) {
    Announce("GAUSSIAN realizations: adding mean values to pixels... ");
    for (i=0; i<Nfields; i++) {
      if (fieldlist.mean(i)!=0.0) {
#pragma omp parallel for 
	for(j=0; j<npixels; j++) mapf[i][j] = mapf[i][j] + fieldlist.mean(i);
      }
    }
    Announce();
  }


  /*** Part 5.2: Generate convergence maps by line of sight (LoS) integration ***/
  
  // If requested, integrate density along the LoS to get convergence:
  if(config.readi("DENS2KAPPA")==1) {
    cout << "Will perform LoS integration over density fields:\n";
    double **KappaWeightTable;
    Healpix_Map<MAP_PRECISION> *IntDens, *tempmapf;
    int zsource, Nintdens=0;

    // Error checking (density fields must have continuous redshift coverage):
    k = fieldlist.CheckZ4Int();
    cout << "   Found "<<k<<" density fields.\n";
    if (k==0) error("flask: no density field found for integrating");
    
    // Compute Kernel:
    Announce("   Tabulating integration kernel... ");
    KappaWeightTable = matrix<double>(0, Nfields-1, 0, Nfields-1);
    for (i=0; i<Nfields; i++) 
      for (j=0; j<Nfields; j++) 
	KappaWeightTable[i][j] = AvgKappaWeightByZ(cosmo, fieldlist.zmin(j), fieldlist.zmax(j), fieldlist.zmax(i)) 
	  * (fieldlist.zmax(j)-fieldlist.zmin(j));
    
    //TabulateKappaWeight(KappaWeightTable, cosmo, fieldlist);
    Announce();
    
    // Do the integration:
    Announce("   Integrating densities... ");
    IntDens = vector<Healpix_Map <MAP_PRECISION> >(0,Nfields-1);
    for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==fgalaxies) { // LOOP over galaxy fields and redshift bins (as sources).
	Nintdens++;
	IntDens[i].SetNside(nside, RING); IntDens[i].fill(0);   // Allocate map at intdens(f,z=z_source) to receive delta(f,z) integral. 
	fieldlist.Index2fFixed(i, &f, &zsource);
#pragma omp parallel for private(z, m)
	for (j=0; j<npixels; j++) {                             // LOOP over pixels (lines of sight).
	  for (z=0; z<=zsource; z++) {                          // LOOP over redshift z (integrating).
	    m = fieldlist.fFixedIndex(f, z);
	    IntDens[i][j] += KappaWeightTable[i][m]*mapf[m][j]; // Sum contributions in the same pixel. 
	  }
	}
      }
    free_matrix(KappaWeightTable, 0, Nfields-1, 0, Nfields-1);
    Announce();

    // Print table with integrated densities statistics:
    filename = config.reads("DENS2KAPPA_STAT");
    if (filename!="0") {
      if (filename=="1") {
	Announce("   Computing integrated density statistics... ");
	cout << endl;
	PrintMapsStats(IntDens, fieldlist, lognormal);
	cout << endl;
	Announce();
      }
      else {
	Announce("   Computing integrated density statistics... ");
	outfile.open(filename.c_str());
	if (!outfile.is_open()) warning("flask: cannot open file "+filename);
	PrintMapsStats(IntDens, fieldlist, lognormal, &outfile);
	outfile.close();
	Announce();
      	cout << ">> DENS2KAPPA_STAT written to "+filename<<endl;
      }
      
    }
    if (ExitAt=="DENS2KAPPA_STAT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
    }

    
    // Join Integrated density to other maps:
    Announce("   Concatenating integrated density data to main data...");
    int *ftemp, *fName, *zName;
    double **ztemp, *mtemp, *stemp;

    // Allocate temporary memory:
    fName    = vector<int>   (0, Nfields+Nintdens-1);
    zName    = vector<int>   (0, Nfields+Nintdens-1);
    ftemp    = vector<int>   (0, Nfields+Nintdens-1);
    mtemp    = vector<double>(0, Nfields+Nintdens-1);
    stemp    = vector<double>(0, Nfields+Nintdens-1);
    ztemp    = matrix<double>(0, Nfields+Nintdens-1, 0, 1);
    tempmapf = vector<Healpix_Map<MAP_PRECISION> >(0, Nfields+Nintdens-1);	  
    // Copy original maps and integrated density maps (and their infos) to same arrays:
    k=0;
    for(i=0; i<Nfields; i++) {
      // Copy original:
      fieldlist.Index2Name(i, &(fName[i]), &(zName[i]));
      ftemp[i]    = fieldlist.ftype(i);
      mtemp[i]    = fieldlist.mean(i);
      stemp[i]    = fieldlist.shift(i);
      ztemp[i][0] = fieldlist.zmin(i);
      ztemp[i][1] = fieldlist.zmax(i);
      tempmapf[i].SetNside(nside, RING);
      tempmapf[i].Import(mapf[i]);
      mapf[i].SetNside(1, RING);
      // Copy integrated densities:
      if (fieldlist.ftype(i)==fgalaxies) {
	k++;
	fName[Nfields-1+k]    = Nf + fName[i];
	zName[Nfields-1+k]    =      zName[i];
	ftemp[Nfields-1+k]    = flensing;            // NOTE: means and shifts are not being set for integrated densities.     
	ztemp[Nfields-1+k][0] = fieldlist.zmax(i);   // The convergence from integration applies to sources
	ztemp[Nfields-1+k][1] = fieldlist.zmax(i);   // located sharply at the end of the bin.
	tempmapf[Nfields-1+k].SetNside(nside, RING);
	tempmapf[Nfields-1+k].Import(IntDens[i]);
	IntDens[i].SetNside(1, RING);
      }   
    }
    // Pass restructured data to main variables:
    free_vector(mapf,    0, Nfields-1); mapf = tempmapf;
    free_vector(IntDens, 0, Nfields-1);
    fieldlist.Build(fName, zName, Nfields+Nintdens, ftemp, ztemp, mtemp, stemp);
    free_vector(fName, 0, Nfields+Nintdens-1);
    free_vector(zName, 0, Nfields+Nintdens-1);
    free_vector(ftemp, 0, Nfields+Nintdens-1);
    free_vector(mtemp, 0, Nfields+Nintdens-1);
    free_vector(stemp, 0, Nfields+Nintdens-1);
    free_matrix(ztemp, 0, Nfields+Nintdens-1, 0, 1);
    Nfields = fieldlist.Nfields(); 
    Nf      = fieldlist.Nfs(); 
    Nz      = fieldlist.Nzs();
    Announce();
  } // End of IF compute convergence by density LoS integration.
  else if (config.readi("DENS2KAPPA")!=0) warning("flask: unknown DENS2KAPPA option: skipping density LoS integration.");
  
  // Write final map to file as a table if requested:
  GeneralOutput(mapf, config, "MAP_OUT", fieldlist);
  // Map output to fits and/or tga files:
  GeneralOutput(mapf, config, "MAPFITS_PREFIX", fieldlist, 1);
  
  // Exit if this is the last output requested:
  if (ExitAt=="MAP_OUT" ||
      ExitAt=="MAPFITS_PREFIX") {
    PrepareEnd(StartAll); return 0;
  }

  // If requested, recover alms and/or Cls from maps:
  RecoverAlmCls(mapf, fieldlist, "RECOVALM_OUT", "RECOVCLS_OUT", config);
  // Exit if this is the last output requested:
  if (ExitAt=="RECOVALM_OUT" || ExitAt=="RECOVCLS_OUT") {
    PrepareEnd(StartAll); return 0;
  }
    

  /*** Part 5.3: Compute shear maps if necessary ***/

  Healpix_Map<MAP_PRECISION> *gamma1f, *gamma2f;    

  if (yesShear==1) {
    gamma1f = vector<Healpix_Map <MAP_PRECISION> >(0,Nfields-1);
    gamma2f = vector<Healpix_Map <MAP_PRECISION> >(0,Nfields-1);
    lmax    = config.readi("SHEAR_LMAX");                            // ATTENTION! Redefining lmax here!
    Alm<xcomplex <ALM_PRECISION> > Eflm(lmax,lmax), Bflm(lmax,lmax); // Temp memory
    arr<double> weight(2*mapf[0].Nside());                           // Temp memory
    for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) {                    // B-modes are zero for weak lensing.
#if USEXCOMPLEX // For compatibility with Healpix versions <=3.20 and >=v3.30.          
	Bflm(l,m).Set(0,0);
#else
	Bflm(l,m).real(0);  
	Bflm(l,m).imag(0);
#endif
      }

    // LOOP over convergence fields:
    for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==flensing) {
	fieldlist.Index2Name(i, &f, &z);
	cout << "** Will compute shear for f"<<f<<"z"<<z<<":\n";
      
	// Preparing memory:
	Announce("   Allocating and cleaning memory... ");
	gamma1f[i].SetNside(nside, RING); gamma1f[i].fill(0);
	gamma2f[i].SetNside(nside, RING); gamma2f[i].fill(0);
	Announce();
      
	// LOGNORMAL REALIZATIONS: get convergence alm's from lognormal convergence map:
	if (dist==lognormal) {
	  PrepRingWeights(1, weight, config);
	  Announce("   Transforming convergence map to harmonic space... ");
	  if (lmax>nside) warning("SHEAR_LMAX > NSIDE introduces noise in the transformation.");
	  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) {
#if USEXCOMPLEX // For compatibility with Healpix versions <=3.20 and >=v3.30.          
	      Eflm(l,m).Set(0,0);
#else
	      Eflm(l,m).real(0);
	      Eflm(l,m).imag(0);
#endif
	    }
	  map2alm_iter(mapf[i], Eflm, 1, weight); // Get klm.
	  Announce();
	}
 
	if (dist!=homogeneous) {
	  // Calculate shear E-mode alm's from convergence alm's:
	  Announce("   Computing shear harmonic coefficients from klm... ");
	  if (dist==lognormal)     Kappa2ShearEmode(Eflm, Eflm);
	  else if (dist==gaussian) Kappa2ShearEmode(Eflm, aflm[i]);
	  Announce();
	}
	else {
	  Announce("HOMOGENEOUS realizations: setting shear E-mode to zero... ");
	  for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) {
#if USEXCOMPLEX // For compatibility with Healpix versions <=3.20 and >=v3.30.          
	      Eflm(l,m).Set(0,0);
#else
	      Eflm(l,m).real(0);
	      Eflm(l,m).imag(0);
#endif
	    }
	  Announce();
	}
	fieldlist.Index2Name(i, &f, &z);
	GeneralOutput(Eflm, config, "SHEAR_ALM_PREFIX", f, z);

	// Go from shear E-mode alm's to gamma1 and gamma2 maps:
	Announce("   Transforming harmonic coefficients into shear map... ");
	alm2map_spin(Eflm, Bflm, gamma1f[i], gamma2f[i], 2);
	Announce();
	// Write kappa, gamma1 and gamma2 to FITS file:
	fieldlist.Index2Name(i, &f, &z);
	GeneralOutput(mapf[i], gamma1f[i], gamma2f[i], config, "SHEAR_FITS_PREFIX", f, z);

      } // End of LOOP over convergence fields.

    // Memory deallocation:
    if (dist==gaussian && yesShear==1) free_vector(aflm, 0, Nfields-1);
    weight.dealloc();
    Eflm.Set(0,0); 
    Bflm.Set(0,0);

    // Exit if this is the last output requested:
    if (ExitAt=="SHEAR_ALM_PREFIX"  ||
	ExitAt=="SHEAR_FITS_PREFIX") {
      PrepareEnd(StartAll); return 0;
    }

    // Output shear maps to TEXT tables:
    GeneralOutput(gamma1f, gamma2f, config, "SHEAR_MAP_OUT", fieldlist);
    // Exit if this is the last output requested:
    if (ExitAt=="SHEAR_MAP_OUT") {
       PrepareEnd(StartAll); return 0;
    }
  } // End of IF we should compute shear.
  


  /************************************/
  /*** Part 6: Maps to Observables  ***/
  /************************************/
  
  SelectionFunction selection;
  MAP_PRECISION maskval;

  // Read in selection functions from FITS files and possibly text files (for radial part):
  Announce("Reading selection functions from files... ");
  selection.load(config, fieldlist); 
  if (selection.Nside()!=-2 && selection.Nside()!=mapf[0].Nside())
    error("flask: Selection function and maps have different number of pixels.");
  if (selection.Scheme()!=-2 && selection.Scheme()!=mapf[0].Scheme()) 
    error("flask: Selection function and maps have different pixel ordering schemes.");
  Announce();
  // Set value used for masked region:
  if (config.readi("USE_UNSEEN")==1) maskval = UNSEEN;
  else if (config.readi("USE_UNSEEN")==0) maskval = 0;
  else {maskval=0;  warning("Unknown option for USE_UNSEEN. Setting it to 0.");}


  /*** Galaxy fields ***/
 
  double dw = 1.4851066049791e8/npixels; // Pixel solid angle in arcmin^2.
  int *counter;
  counter = vector<int>(1, MaxThreads);

  // Poisson Sampling the galaxy fields:
  if (config.readi("POISSON")==1) {
    // LOOP over galaxy fields:
    for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==fgalaxies) {
	fieldlist.Index2Name(i, &f, &z);
	sprintf(message, "Poisson sampling f%dz%d... ", f, z); filename.assign(message); 
	Announce(filename);
	for(k=1; k<=MaxThreads; k++) counter[k]=0;
	// LOOP over pixels of field 'i':
#pragma omp parallel for schedule(static) private(k)
	for(j=0; j<npixels; j++) {
	  k = omp_get_thread_num()+1;
	  if (mapf[i][j] < -1.0) { counter[k]++; mapf[i][j]=0.0; } // If density is negative, set it to zero.
	  if (selection.MaskBit(i,j)==1 || selection.MaskBit(i,j)==3 || selection.MaskBit(i,j)==2) mapf[i][j]=maskval;  
	  else mapf[i][j] = gsl_ran_poisson(rnd[k], selection(i,j)*(1.0+mapf[i][j])*dw);
	}
	Announce();
	j=0; for (k=1; k<=MaxThreads; k++) j+=counter[k];
	cout << "Negative density fraction (that was set to 0): "<<std::setprecision(2)<< ((double)j)/npixels*100 <<"%\n";
      }
  }

  // OR Gaussian Sampling the galaxy fields:
  else if (config.readi("POISSON")==2) {
    // LOOP over galaxy fields:
    for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==fgalaxies) {
	fieldlist.Index2Name(i, &f, &z);
	sprintf(message, "Gaussian sampling f%dz%d... ", f, z); filename.assign(message); 
	Announce(filename);
	for(k=1; k<=MaxThreads; k++) counter[k]=0;
	// LOOP over pixels of field 'i':
#pragma omp parallel for schedule(static) private(k)
	for(j=0; j<npixels; j++) {
	  k = omp_get_thread_num()+1;
	  if (selection.MaskBit(i,j)==1 || selection.MaskBit(i,j)==3 || selection.MaskBit(i,j)==2) mapf[i][j]=maskval;
	  else {
	    mapf[i][j] = selection(i,j)*(1.0+mapf[i][j])*dw + gsl_ran_gaussian(rnd[k], sqrt(selection(i,j)*dw)); 
	    if (mapf[i][j] < 0) counter[k]++; // Count pixels that have negative number of galaxies after sampling.
	  }
	}
	Announce();
	j=0; for (k=1; k<=MaxThreads; k++) j+=counter[k];
	cout << "Negative counts fraction (after sampling): "<<std::setprecision(2)<< ((double)j)/npixels*100 <<"%\n";
      }
  }
  
  // Just generate the expected number density, if requested:
  else if (config.readi("POISSON")==0) {
    // LOOP over galaxy fields:
    for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==fgalaxies) {
	fieldlist.Index2Name(i, &f, &z);
	sprintf(message,"Using expected (no sampling) number counts for f%dz%d...", f, z); filename.assign(message);
	Announce(message);	
	for(k=1; k<=MaxThreads; k++) counter[k]=0;
	// LOOP over pixels of field 'i':
#pragma omp parallel for schedule(static) private(k)
	for(j=0; j<npixels; j++) {
	  k = omp_get_thread_num()+1;
	  if (selection.MaskBit(i,j)==1 || selection.MaskBit(i,j)==3 || selection.MaskBit(i,j)==2) mapf[i][j]=maskval;
	  else {
	    mapf[i][j] = selection(i,j)*(1.0+mapf[i][j])*dw;
	    if (mapf[i][j] < 0) counter[k]++; // Count pixels that have negative number of galaxies.
	  }
	}
	Announce();
	j=0; for (k=1; k<=MaxThreads; k++) j+=counter[k];
	cout << "Negative counts fraction: "<<std::setprecision(2)<< ((double)j)/npixels*100 <<"%\n";
      }
  }
  
  else error ("flask: unknown POISSON option.");
  free_vector(counter, 1, MaxThreads);
  

  /*** Convergence fields ***/

  // LOOP over fields:
  if (config.reads("MAPWERFITS_PREFIX")!="0" || config.reads("MAPWER_OUT")!="0" || config.reads("ELLIPFITS_PREFIX")!="0") {
    if (CountGalaxyFields(fieldlist)==0 || config.readi("SELEC_SEPARABLE")==1) {
      for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==flensing) {
	  fieldlist.Index2Name(i, &f, &z);
	  sprintf(message, "Masking f%dz%d... ", f, z); filename.assign(message); 
	  Announce(filename);
#pragma omp parallel for
	  for(j=0; j<npixels; j++)
	    // No noise is added, we only apply the mask:
	    if (selection.MaskBit(i,j)==1 || selection.MaskBit(i,j)==3 || selection.MaskBit(i,j)==2) mapf[i][j]=maskval;
	  Announce();
	}
    }
  }
  
  /*** Output of observable convergence and galaxy count maps ***/
  
  // Write final map to file as a table if requested:
  GeneralOutput(mapf, config, "MAPWER_OUT", fieldlist);
  // Map output to fits and/or tga files:
  GeneralOutput(mapf, config, "MAPWERFITS_PREFIX", fieldlist, 1);
  
  // Exit if this is the last output requested:
  if (ExitAt=="MAPWER_OUT" ||
      ExitAt=="MAPWERFITS_PREFIX") {
    PrepareEnd(StartAll); return 0;
  }

  
  /*** Ellipticity fields ***/
  
  if (config.reads("ELLIP_MAP_OUT")!="0" || config.reads("ELLIPFITS_PREFIX")!="0") {
    Healpix_Map <MAP_PRECISION> *e1Mapf, *e2Mapf;
    e1Mapf = vector<Healpix_Map <MAP_PRECISION> >(0,Nfields-1);
    e2Mapf = vector<Healpix_Map <MAP_PRECISION> >(0,Nfields-1);
    esig   = config.readd("ELLIP_SIGMA");
    i = CountGalaxyFields(fieldlist);
    
    // If there are galaxy fields, use them to create noise for the ellipticity maps:
    if (i>0) {
      cout << "Will compute ellip. for z slices w/ gals. and lens. fields:\n";
      // Find out which is the lensing field to use:
      k = CountLensingFields(fieldlist);
      if (k>1) warning("flask: found multiple lensing fields, will use the last one found.");
      for(f=0, m=0; f<fieldlist.Nfs() && m<k; f++) if (fieldlist.ftype(fieldlist.fFixedIndex(f, 0))==flensing) m++;
      fieldlist.fFixedName(f-1,0, &k, &l);
      // LOOP over galaxy fields:
      for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==fgalaxies) {
	  fieldlist.Index2Name(i, &f, &z);
	  fieldlist.Name2Index(k, z, &j, 0);  // << Get index of the corresponding lensing field.
	  if (j>-1) {                         // << If the lensing field exists, get the corresponding ellipticities:
	    e1Mapf[i].SetNside(nside,RING);  e2Mapf[i].SetNside(nside,RING);
	    sprintf(message, "Generating ellipticity for f%dz%d...", f, z); filename.assign(message); 
	    Announce(filename);
#pragma omp parallel for schedule(static) private(k)
	    for(m=0; m<npixels; m++) {
	      k = omp_get_thread_num()+1;
	      // Mask ellipticity map at galaxy-free pixels:
	      if (mapf[i][m]<=0) { e1Mapf[i][m]=maskval;  e2Mapf[i][m]=maskval; }
	      // If galaxies are present, average ellipticity in pixel is the reduced shear plus the average of the error:
	      else {
		if (esig>0.0) {
		  e1Mapf[i][m] = gamma1f[j][m]/(1.0-mapf[j][m]) + gsl_ran_gaussian(rnd[k], esig/sqrt(mapf[i][m]));
		  e2Mapf[i][m] = gamma2f[j][m]/(1.0-mapf[j][m]) + gsl_ran_gaussian(rnd[k], esig/sqrt(mapf[i][m]));
		}
		else {
		  e1Mapf[i][m] = gamma1f[j][m]/(1.0-mapf[j][m]);
		  e2Mapf[i][m] = gamma2f[j][m]/(1.0-mapf[j][m]);
		}
	      }
	    }
	    Announce();
	    // Write counts, ellip1 and ellip2 to FITS file:
	    GeneralOutput(mapf[j], e1Mapf[i], e2Mapf[i], config, "ELLIPFITS_PREFIX", f, z);
	  }
	} // End of LOOP over galaxy fields.
    }
    
    // If not, use the shear sources selection function that should have been provided:
    else {
      cout << "Will compute ellip. from lens. fields' selection functions:\n";
      // LOOP over lensing fields:
      for (i=0; i<Nfields; i++) if (fieldlist.ftype(i)==flensing) {
	  fieldlist.Index2Name(i, &f, &z);
	  e1Mapf[i].SetNside(nside,RING);  e2Mapf[i].SetNside(nside,RING);
	  sprintf(message, "Generating ellipticity for f%dz%d...", f, z); filename.assign(message); 
	  Announce(filename);
#pragma omp parallel for schedule(static) private(k)
	  for(m=0; m<npixels; m++) {
	    k = omp_get_thread_num()+1;
	    // Mask ellipticity map at masked regions:
	    if (selection(i,m)<=0) { e1Mapf[i][m]=maskval;  e2Mapf[i][m]=maskval; }
	    // In unmasked regions, average ellipticity in pixel is the reduced shear plus the average of the error:
	    else {
	      if (esig>0.0) {
		e1Mapf[i][m] = gamma1f[i][m]/(1.0-mapf[i][m]) + gsl_ran_gaussian(rnd[k], esig/sqrt(selection(i,m)*dw));
		e2Mapf[i][m] = gamma2f[i][m]/(1.0-mapf[i][m]) + gsl_ran_gaussian(rnd[k], esig/sqrt(selection(i,m)*dw));
	      }
	      else {
		e1Mapf[i][m] = gamma1f[i][m]/(1.0-mapf[i][m]);
		e2Mapf[i][m] = gamma2f[i][m]/(1.0-mapf[i][m]);
	      }
	    }
	  }
	  Announce();
	  // Write counts, ellip1 and ellip2 to FITS file:
	  GeneralOutput(mapf[i], e1Mapf[i], e2Mapf[i], config, "ELLIPFITS_PREFIX", f, z);
	} // End of LOOP over lensing fields.
    }

    // Output ellipticity maps to TEXT tables:
    GeneralOutput(e1Mapf, e2Mapf, config, "ELLIP_MAP_OUT", fieldlist);

    free_vector(e1Mapf, 0, Nfields-1);
    free_vector(e2Mapf, 0, Nfields-1);
  } // End of compute ellipticity maps.
  
  // Exit if this is the last output requested:
  if (ExitAt=="ELLIP_MAP_OUT" ||
      ExitAt=="ELLIPFITS_PREFIX") {
    PrepareEnd(StartAll); return 0;
  }


  /**********************************/
  /*** Part 7: Generate catalog   ***/
  /**********************************/

  CAT_PRECISION **catalog;
  char **catSet;
  double ellip1, ellip2, randz, rdist;
  int gali, cellNgal, ncols, AngularCoord;
  long *ThreadNgals, Ngalaxies, kl, Ncells, PartialNgal, longNz;  
  pointing ang;
  int ziter, fiter;
  std::string CatalogHeader;
  int theta_pos, phi_pos, z_pos, r_pos, galtype_pos, kappa_pos, gamma1_pos, gamma2_pos, 
    ellip1_pos, ellip2_pos, pixel_pos, maskbit_pos, ra_pos, dec_pos;

  esig = config.readd("ELLIP_SIGMA");
  
  Announce("Counting galaxies... ");
  // Partial count of galaxies for each thread (assuming the threads get the data in order):
  ThreadNgals = vector<long>(0, MaxThreads);       
  for(l=0; l<=MaxThreads; l++) ThreadNgals[l] = 0;
  longNz  = (long)Nz;
  Ncells  = ((long)npixels)*longNz; // Number of 3D cells.
  // Loop over 3D cells:
#pragma omp parallel for schedule(static) private(l, j, ziter, fiter, i)
  for (kl=0; kl<Ncells; kl++) {
    l     = omp_get_thread_num();
    j     = kl/longNz;     // angular pixel.
    ziter = kl%longNz;     // z slice.
    for (fiter=0; fiter<fieldlist.Nf4z(ziter); fiter++) {
      i = fieldlist.zFixedIndex(fiter, ziter);
      if (fieldlist.ftype(i)==fgalaxies && mapf[i][j]>0) ThreadNgals[l+1] += (long)mapf[i][j];
    }
  }
  // Compute the cummulative sum by thread, which is the catalogue initial position for thar thread:
  for (l=2; l<=MaxThreads; l++) ThreadNgals[l] += ThreadNgals[l-1];
  Ngalaxies = ThreadNgals[MaxThreads];
  Announce();     
  cout << "# of galaxies: "<<Ngalaxies<<endl;
  if (Ngalaxies>INT_MAX) warning("flask: catalogue generation not tested for this amount of galaxies");
  
  // Using transposed catalog (catalog[col][row]), better for FITS outputting:
  CatalogHeader = config.reads("CATALOG_COLS");
  ncols         = CountWords(CatalogHeader);
  catalog       = matrix<CAT_PRECISION>(0,ncols-1, 0,Ngalaxies-1); 
  catSet        = matrix<char>         (0,ncols-1, 0,Ngalaxies-1);
  for (kl=0; kl<Ngalaxies; kl++) for (j=0; j<ncols; j++) catSet[j][kl]=0;
  
  // Find position of entries according to catalog header:
  theta_pos    = GetSubstrPos("theta"  , CatalogHeader); 
  phi_pos      = GetSubstrPos("phi"    , CatalogHeader);
  ra_pos       = GetSubstrPos("ra"     , CatalogHeader); 
  dec_pos      = GetSubstrPos("dec"    , CatalogHeader);
  z_pos        = GetSubstrPos("z"      , CatalogHeader);  
  r_pos        = GetSubstrPos("r"      , CatalogHeader);  
  galtype_pos  = GetSubstrPos("galtype", CatalogHeader);  
  kappa_pos    = GetSubstrPos("kappa"  , CatalogHeader);  
  gamma1_pos   = GetSubstrPos("gamma1" , CatalogHeader);  
  gamma2_pos   = GetSubstrPos("gamma2" , CatalogHeader);  
  ellip1_pos   = GetSubstrPos("ellip1" , CatalogHeader);  
  ellip2_pos   = GetSubstrPos("ellip2" , CatalogHeader);  
  pixel_pos    = GetSubstrPos("pixel"  , CatalogHeader); 
  maskbit_pos  = GetSubstrPos("maskbit", CatalogHeader);

  // Allow Change of Coordinates if RA and DEC were set as catalog columns:
  // For the catalog, ra, dec, theta, phi in CATALOG_COLS overrides ANGULAR_COORD. 
  AngularCoord = config.readi("ANGULAR_COORD");
  OrganizeAngularCoord(&AngularCoord, &phi_pos, &theta_pos, &ra_pos, &dec_pos, CatalogHeader);

  // Warning against multiple or none lensing fields at the same redshift:
  //k=0;
  //for(f=0; f<Nf; f++) if (fieldlist.ftype(fieldlist.fFixedIndex(f, 0))==flensing) k++;
  k = CountLensingFields(fieldlist);
  if (k>1) warning("flask: found multiple lensing fields, will leave last one in catalogue.");
  if (k<1 && (kappa_pos!=-1 || gamma1_pos!=-1 || gamma2_pos!=-1 || ellip1_pos!=-1 || ellip1_pos!=-1)) 
    warning("flask: missing lensing information required to build catalogue.");

  // LOOP over 3D cells (pixels and redshifts):
  Announce("Generating catalog... ");
  if (r_pos!=-1) ComDist(cosmo,1.0);   // Initialize Comoving distance formula if requested.
  PartialNgal=0;                       // Counter of galaxies inside thread.
#pragma omp parallel for schedule(static) private(l, j, ziter, gali, fiter, i, m, ang, ellip1, ellip2, randz, rdist, cellNgal, f, z) firstprivate(PartialNgal)
  // Since this FOR has the same parameters as the one above for counting, thread assignment should be the same. 
  for (kl=0; kl<Ncells; kl++) {
    l        = omp_get_thread_num();   // Processor number.
    j        = kl/longNz;              // pixel.
    ziter    = kl%longNz;              // z slice.
    gali     = 0;                      // Galaxy number inside cell.
    cellNgal = 0;                      // Total galaxy number inside cell.
    
    // Count total number of galaxies of all types inside cell:
    for (fiter=0; fiter<fieldlist.Nf4z(ziter); fiter++) {
      i = fieldlist.zFixedIndex(fiter, ziter);
      if (fieldlist.ftype(i)==fgalaxies && mapf[i][j]>0) cellNgal += (int)mapf[i][j];
    }
 
    // LOOP over field IDs:
    for (fiter=0; fiter<fieldlist.Nf4z(ziter); fiter++) {
      i = fieldlist.zFixedIndex(fiter, ziter);
      fieldlist.Index2Name(i, &f, &z);
      
      // Add entry of type GALAXY:      
      if (fieldlist.ftype(i)==fgalaxies) for(m=0; m<(int)mapf[i][j]; m++) {
	  if (theta_pos!=-1 || phi_pos!=-1) ang   = RandAngInPix(rnd[l+1], mapf[i], j);
	  if (z_pos!=-1 || r_pos!=-1)       randz = selection.RandRedshift(rnd[l+1],i,j);
	  if (maskbit_pos!=-1)              k     = selection.MaskBit(i,j);
	  if (r_pos!=-1)                    rdist = ComDist(cosmo, randz);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, theta_pos  , ang.theta, catSet);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, phi_pos    , ang.phi  , catSet);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, z_pos      , randz    , catSet);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, r_pos      , rdist    , catSet);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, galtype_pos, f        , catSet);	    
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, pixel_pos  , j        , catSet);
	  CatalogFill(catalog, ThreadNgals[l]+PartialNgal+gali, maskbit_pos, k        , catSet);
	  gali++;
	}
      
      // Add entry of type SHEAR:
      else if (fieldlist.ftype(i)==flensing) for (m=0; m<cellNgal; m++) {
	  CatalogFill  (catalog, ThreadNgals[l]+PartialNgal+m, kappa_pos , mapf[i][j]   , catSet);
	  if (yesShear==1) {
	    if (ellip1_pos!=-1 || ellip2_pos!=-1) GenEllip(rnd[l+1], esig, mapf[i][j], gamma1f[i][j], gamma2f[i][j], &ellip1, &ellip2);
	    CatalogFill(catalog, ThreadNgals[l]+PartialNgal+m, gamma1_pos, gamma1f[i][j], catSet);
	    CatalogFill(catalog, ThreadNgals[l]+PartialNgal+m, gamma2_pos, gamma2f[i][j], catSet);
	    CatalogFill(catalog, ThreadNgals[l]+PartialNgal+m, ellip1_pos, ellip1       , catSet);
	    CatalogFill(catalog, ThreadNgals[l]+PartialNgal+m, ellip2_pos, ellip2       , catSet);  
	  }
	}
    } // End of LOOP over field IDs.
    PartialNgal += (long)cellNgal;
  } // End of LOOP over cells.  
  free_vector(ThreadNgals, 0, MaxThreads);

  // Check if every entry was set once and only once:
  for (kl=0; kl<Ngalaxies; kl++) for (j=0; j<ncols; j++) {
      if (catSet[j][kl]<1) {
	printf("j=%d kl=%ld N=%d\n", j, kl, catSet[j][kl]);
	warning("flask: Catalog has missing information.");}
      if (catSet[j][kl]>1) {
	printf("j=%d kl=%ld N=%d\n", j, kl, catSet[j][kl]);
	warning("flask: Catalog entry being set more than once.");}
    }
  free_matrix(catSet, 0,ncols-1, 0,Ngalaxies-1);
  Announce();

  // Memory deallocation: all required info from now on should be in catalog.
  Announce("Deallocating maps and fields info... ");
  free_vector(mapf,    0, Nfields-1 );
  if (yesShear==1) free_vector(gamma1f, 0, Nfields-1);
  if (yesShear==1) free_vector(gamma2f, 0, Nfields-1);
  for (i=0; i<=MaxThreads; i++) gsl_rng_free(rnd[i]);
  free_vector(rnd, 0,MaxThreads+1);
  Announce();

  // Change angular coordinates if requested:
  ChangeCoord(catalog, theta_pos, phi_pos, Ngalaxies, AngularCoord);

  // Write catalog to file if requested:
  if (config.reads("CATALOG_OUT")!="0") {
    filename = config.reads("CATALOG_OUT");
    k        = FileFormat(filename);
    switch (k) {
      // TEXT file:
    case ascii_format: 
      outfile.open(filename.c_str());
      if (!outfile.is_open()) warning("flask: cannot open file "+filename);
      else {
	outfile << "# "<< CatalogHeader <<endl;
	PrintVecs(catalog, Ngalaxies, ncols, &outfile); 
	outfile.close();
	cout << ">> Catalog written to " << filename << endl;
      }
      break;
      // FITS file:
    case fits_format: 
      WriteCatalog2Fits(filename, catalog, Ngalaxies, config);
      cout << ">> Catalog written to " << filename << endl;
      break;
      // Unknown: 
    case unknown_format: 
      warning("flask: unknown catalogue file format, no output performed");
      break;
      // Weird: 
    default:
      warning("flask: uninplemented catalogue file format, check code");
      break;
    }
  }  
  free_matrix(catalog, 0,ncols-1, 0,Ngalaxies-1);
  

  // End of the program
  PrepareEnd(StartAll); return 0;
}
