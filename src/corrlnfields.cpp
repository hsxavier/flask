/* corrlnfields: Written by Henrique S. Xavier on Nov-2014
   e-mail: hsxavier@if.usp.br
 */

#include <iostream>
#include "corrlnfields_aux.hpp" // Auxiliary functions made for this program.
#include "GeneralOutput.hpp"    // Various file output functions.
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include "s2kit10_naive.hpp"    // For Discrete Legendre Transforms.
#include "Cosmology.hpp"        // Parameters and formulas.
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
#include "SelectionFunc.hpp"
#include "RegularizeCov.hpp"
#include <gsl/gsl_eigen.h> // debug

/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl;                     // Basic stuff.
  using namespace definitions;                          // Global definitions.
  using namespace ParDef; ParameterList config;         // Easy configuration file use.
  Cosmology cosmo;                                      // Cosmological parameters.
  char message[100];                                    // Handling warnings and errors.
  std::string filename, tempstr;
  std::ofstream outfile;                                // File for output.
  enum simtype {gaussian, lognormal}; simtype dist;     // For specifying simulation type.
  gsl_matrix **CovByl; 
  int status, i, j, l, m, N1, N2, Nfields, mmax, *ftype;
  double *means, *shifts, **zrange; 
  long long1, long2;
  FILE* stream; int NinputCls; std::string *filelist;                       // To list input Cls.
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


  /********************************************/
  /*** PART 2: Load C(l)s and organize them ***/
  /********************************************/
  const int HWMAXL = 10000000; int lastl = HWMAXL, Nls;
  int a1, a2, b1, b2, Nlinput, **fnz, **NentMat;
  long *Nentries, ncols;
  double ***ll, ***Cov, *wrapper[2];
  bool **IsSet;
  
  // Listing files to use based on CL_PREFIX:
  // Find out how many input C(l)'s we have.
  cout << "Will load input C(l)s:\n";
  sprintf(message, "ls %s* | wc -l", config.reads("CL_PREFIX").c_str()); 
  stream = popen(message, "r");
  if ((stream=popen(message, "r")) == NULL) error("corrlnfields: cannot 'ls | wc' output");
  fscanf(stream, "%d", &NinputCls); pclose(stream);
  // Get list of input C(l) files:
  filelist = vector<std::string>(0,NinputCls-1);
  sprintf(message, "ls %s*", config.reads("CL_PREFIX").c_str());
  if ((stream=popen(message, "r")) == NULL) error("corrlnfields: cannot pipe 'ls' output");
  for (i=0; i<NinputCls; i++) {
    fscanf(stream, "%s", message);
    filelist[i].assign(message);
  }
  pclose(stream);

  // Get file list and find out how many C(l)s there are:  
  N1=0; N2=0; Nlinput=0;
  Nentries = vector<long>(0,NinputCls-1);
  for (i=0; i<NinputCls; i++) {
    getcovid(filelist[i], &a1, &a2, &b1, &b2);
    if (a1>N1) N1=a1; if (b1>N1) N1=b1;                // Get number of fields.
    if (a2>N2) N2=a2; if (b2>N2) N2=b2;                // Get number of z bins.
    CountEntries(filelist[i], &(Nentries[i]), &ncols); // Get number of Nls.
    if (ncols!=2) error("corrlnfields: wrong number of columns in file "+filename);
    if (Nentries[i]>Nlinput) Nlinput=Nentries[i];          // Record maximum number of ls.
  }
  //Check if number of fields in INFO is the same as in the Cls:
  if (Nfields != N1*N2) error("corrlnfields: number of means and shifts do not match number of C(l)s.");
  

  // Allocate memory to store C(l)s:
  // First two indexes are CovMatrix indexes and last is for ll.
  // fnz stores the order that the fields are stored in CovMatrix.
  fnz     =     matrix<int>(0, Nfields-1, 0, 1);                 // Records what field is stored in each element of CovMatrix.
  fnzSet  =    vector<bool>(0, Nfields-1);                       // For bookkeeping.
  ll      = tensor3<double>(0, Nfields-1, 0, Nfields-1, 0, Nlinput); // Records the ll for each C(l) file. 
  Cov     = tensor3<double>(0, Nfields-1, 0, Nfields-1, 0, Nlinput); // Records the C(l) for each C(l) file.
  IsSet   =    matrix<bool>(0, Nfields-1, 0, Nfields-1);           // For bookkeeping.
  NentMat =     matrix<int>(0, Nfields-1, 0, Nfields-1);           // Number of C(l) entries in file.
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) IsSet[i][j]=0;
  for(i=0; i<Nfields; i++) fnzSet[i]=0;
  
  // Read C(l)s and store in data-cube:
  for (m=0; m<NinputCls; m++) {
    // Find CovMatrix indexes of C(l):
    getcovid(filelist[m], &a1, &a2, &b1, &b2);
    fz2n(a1, a2, &i, N1, N2); fz2n(b1, b2, &j, N1, N2); 
    cout << filelist[m] << " goes to ["<<i<<", "<<j<<"]" << endl;
    // Record the order of the fields in CovMatrix:
    if (fnzSet[i]==0) { fnz[i][0] = a1; fnz[i][1] = a2; fnzSet[i] = 1; }
    else if (fnz[i][0] != a1 || fnz[i][1] != a2) error("corrlnfields: field order in CovMatrix is messed up!"); 
    if (fnzSet[j]==0) { fnz[j][0] = b1; fnz[j][1] = b2; fnzSet[j] = 1; }
    else if (fnz[j][0] != b1 || fnz[j][1] != b2) error("corrlnfields: field order in CovMatrix is messed up!");
    // Import data:
    wrapper[0] = &(ll[i][j][0]);
    wrapper[1] = &(Cov[i][j][0]);
    ImportVecs(wrapper, Nentries[m], 2, filelist[m].c_str());
    NentMat[i][j] = Nentries[m];
    IsSet[i][j]=1; 
  };
  free_vector(Nentries, 0, NinputCls-1);
  free_vector(filelist, 0, NinputCls-1);

  // Check if every field was assigned a position in the CovMatrix:
  for (i=0; i<Nfields; i++) if (fnzSet[i]==0) error("corrlnfields: some position in CovMatrix is unclaimed.");
  free_vector(fnzSet, 0, Nfields-1);
  // If positions are OK and output required, print them out:
  if (config.reads("FLIST_OUT")!="0") {
    outfile.open(config.reads("FLIST_OUT").c_str());
    if (!outfile.is_open()) error("corrlnfields: cannot open FLIST_OUT file.");
    PrintTable(fnz, Nfields, 2, &outfile);
    outfile.close();
    cout << ">> Written field list to "+config.reads("FLIST_OUT")<<endl;
  }
  free_matrix(fnz, 0, Nfields-1, 0,1);
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="FLIST_OUT") {
      cout << "\nTotal number of warnings: " << warning("count") << endl;
      cout<<endl;
      return 0;
  }

  // Look for the maximum l value described by all C(l)s:
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]==1) {
	if (ll[i][j][NentMat[i][j]-1]>HWMAXL) error ("corrlnfields: too high l in C(l)s: increase HWMAXL.");
	if (ll[i][j][NentMat[i][j]-1]<lastl) lastl = (int)ll[i][j][NentMat[i][j]-1];
      }
  Nls=lastl+1; // l=0 is needed for DLT. Nls is known as 'bandwidth' (bw) in s2kit 1.0 code.
  
    
  /*****************************************************************/
  /*** PART 3: Compute auxiliary gaussian C(l)s if LOGNORMAL     ***/
  /*****************************************************************/
  double *tempCl, *LegendreP, *workspace, *xi, lsup, supindex, *theta, *DLTweights, *lls;

  // Allocate gsl_matrices that will receive covariance matrices for each l.
  cout << "Allocating data-cube necessary for Cholesky decomposition... "; cout.flush();
  tempCl = vector<double>(0, lastl);
  CovByl = GSLMatrixArray(Nls, Nfields, Nfields);
  cout << "done.\n";

  if (dist==lognormal) {
    cout << "LOGNORMAL realizations: will compute auxiliary gaussian C(l)s:\n";
    // Loads necessary memory:
    cout << "Allocating extra memory... "; cout.flush();
    workspace  = vector<double>(0, 16*Nls-1);
    LegendreP  = vector<double>(0, 2*Nls*Nls-1);
    xi         = vector<double>(0, 2*Nls-1);
    lls        = vector<double>(0, lastl);
    DLTweights = vector<double>(0, 4*Nls-1);
    // Initialize vectors:
    for (i=0; i<=lastl; i++) lls[i]=(double)i;
    // angle theta is only necessary for output:
    if (config.reads("XIOUT_PREFIX")!="0" || config.reads("GXIOUT_PREFIX")!="0") {
      cout << "Generating table of sampling angles... "; cout.flush();
      theta    = vector<double>(0, 2*Nls-1);
      ArcCosEvalPts(2*Nls, theta);
      for (i=0; i<2*Nls; i++) theta[i] = theta[i]*180.0/M_PI;
      cout << "done.\n";
    } 
    cout << "done.\n";
    // Loads C(l) exponential suppression:
    lsup     = config.readd("SUPPRESS_L");
    supindex = config.readd("SUP_INDEX"); 
    // Load s2kit 1.0 Legendre Polynomials:
    cout << "Generating table of Legendre polynomials... "; cout.flush();
    PmlTableGen(Nls, 0, LegendreP, workspace);
    cout << "done.\n";
    // Compute s2kit 1.0 Discrete Legendre Transform weights:
    cout << "Calculating forward DLT weights... "; cout.flush();
    makeweights(Nls, DLTweights);
    cout << "done.\n";
  }

  // LOOP over all C(l)s already set.
  for(i=0; i<Nfields; i++)
    for(j=0; j<Nfields; j++) 
      if (IsSet[i][j]==1) {
	cout << "** Transforming C(l) in ["<<i<<", "<<j<<"]:\n";
	// Interpolate C(l) for every l; input C(l) might not be like that:
	cout << "   Interpolating input C(l) for all l's... "; cout.flush();
	GetAllLs(ll[i][j], Cov[i][j], NentMat[i][j], tempCl, lastl, config.readi("EXTRAP_DIPOLE"));
	cout << "              done.\n";
	
	if (dist==lognormal) {              /** LOGNORMAL ONLY **/
	  // Compute correlation function Xi(theta):
	  cout << "   DLT (inverse) to obtain the correlation function... "; cout.flush();
	  ModCl4DLT(tempCl, lastl, lsup, supindex);
	  Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	  cout << "  done.\n";
	  if (config.reads("XIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("XIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    cout << ">> Correlation function written to "+filename<<endl;
	  }
	  // Transform Xi(theta) to auxiliary gaussian Xi(theta):
	  cout << "   Computing associated gaussian correlation function... "; cout.flush(); 
	  status=GetGaussCorr(xi, xi, 2*Nls, means[i], shifts[i], means[j], shifts[j]);
	  cout << "done.\n";
	  if (status==EDOM) error("corrlnfields: GetGaussCorr found bad log arguments.");
	  if (i==j && xi[0]<0) warning("corrlnfields: auxiliary field variance is negative.");
	  if (config.reads("GXIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GXIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    cout << ">> Associated Gaussian correlation function written to "+filename<<endl;
	  }
	  // Transform Xi(theta) back to C(l):
	  cout << "   DLT (forward) to obtain the angular power spectrum... "; cout.flush(); 
	  Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	  ApplyClFactors(tempCl, Nls);
	  cout << "done.\n";
	  if (config.reads("GCLOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GCLOUT_PREFIX"), i, j, N1, N2, lls, tempCl, Nls);
	    cout << ">> C(l) for auxiliary Gaussian variables written to "+filename<<endl;
	  }	  
	}                                 /** END OF LOGNORMAL ONLY **/ 
	
	// Save auxiliary C(l):
	for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j]=tempCl[l];		
      } // End of LOOP over C(l)[i,j] that were set.

  // Memory deallocation:
  free_tensor3(Cov,    0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_tensor3(ll,     0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_matrix(NentMat, 0, Nfields-1, 0, Nfields-1);
  //if (config.reads("XIOUT_PREFIX")!="0" || config.reads("GXIOUT_PREFIX")!="0") free_vector(theta, 0, 2*Nls-1);

  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="XIOUT_PREFIX"  || 
      config.reads("EXIT_AT")=="GXIOUT_PREFIX" || 
      config.reads("EXIT_AT")=="GCLOUT_PREFIX") {
    cout << "\nTotal number of warnings: " << warning("count") << endl;
    cout<<endl;
    return 0;
  }
  
  // Set Cov(l)[i,j] = Cov(l)[j,i]
  cout << "Set remaining covariance matrices elements based on symmetry... "; cout.flush(); 
  for(i=0; i<Nfields; i++)
    for(j=0; j<Nfields; j++) 
      if (IsSet[i][j]==0) {
	if (IsSet[j][i]==0) {
	  sprintf(message,"corrlnfields: [%d,%d] could not be set because [%d,%d] was not set.",i,j,j,i);
	  error(message);
	}
	for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j] = CovByl[l]->data[j*Nfields+i];
	IsSet[i][j] = 1;
      }
  cout << "done.\n";
  free_matrix(IsSet, 0, Nfields-1, 0, Nfields-1);
  
  // Output covariance matrices for each l if requested:
  GeneralOutput(CovByl, config, "COVL_PREFIX");
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="COVL_PREFIX") {
    cout << "\nTotal number of warnings: " << warning("count") << endl;
    cout<<endl;
    return 0;
  }


  /***********************************************************/
  /*** PART 3.5: Obtain regularized input Cls if requested ***/
  /***********************************************************/

  if (config.reads("REG_CL_PREFIX")!="0") {
    cout << "** Will compute regularized lognormal Cls:\n";
    if (dist!=lognormal) warning("corrlnfields: option DIST is not LOGNORMAL. Crashing in 3, 2, 1...");

    // Regularizing auxiliary gaussian covariance matrices:
    for (l=1; l<Nls; l++) {
      cout << "   l: "<<l<<endl;
      RegularizeCov(CovByl[l], config);
    }
    
    // LOOP over fields:
    for(i=0; i<Nfields; i++)
      for(j=i; j<Nfields; j++) {
	cout << "** Transforming C(l) in ["<<i<<", "<<j<<"]:\n";
	// Copy the Cl to a vector:
	for (l=0; l<Nls; l++) tempCl[l] = CovByl[l]->data[i*Nfields+j]; // tudo certo.
	// Compute correlation function Xi(theta):
	cout << "   DLT (inverse) to obtain the correlation function...     "; cout.flush();
	ModCl4DLT(tempCl, lastl, -1, -1); // Suppression not needed (it was already suppressed).
	Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	cout << "done.\n";
	// Get Xi(theta) for lognormal variables:
	cout << "   Getting correlation function for lognormal variables... "; cout.flush();
	GetLNCorr(xi, xi, 2*Nls, means[i], shifts[i], means[j], shifts[j]);
	cout << "done.\n";
	// Compute the Cls:
	cout << "   DLT (forward) to obtain the angular power spectrum...   "; cout.flush(); 
	Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	ApplyClFactors(tempCl, Nls, lsup, supindex);
	cout << "done.\n";
	// Output:
	filename=PrintOut(config.reads("REG_CL_PREFIX"), i, j, N1, N2, lls, tempCl, Nls);
	cout << ">> Regularized lognormal C(l) written to "+filename<<endl;
      } 

  } // End of computing regularized lognormal Cls.
  
  // Freeing memory: from now on we only need CovByl, means, shifts, fnz.
  free_vector(tempCl, 0, lastl);
  if (dist==lognormal) {
    cout << "   DLT memory deallocation... "; cout.flush();
    free_vector(workspace, 0, 16*Nls-1);
    free_vector(LegendreP, 0, 2*Nls*Nls-1);
    free_vector(xi, 0, 2*Nls-1);
    free_vector(lls, 0, lastl);
    free_vector(DLTweights, 0, 4*Nls-1); 
    cout << "done.\n";
  }
  
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="REG_CL_PREFIX") {
    cout << "\nTotal number of warnings: " << warning("count") << endl;
    cout<<endl;
    return 0;
  }

  /*********************************************************/
  /*** PART 4: Cholesky decomposition and alm generation ***/
  /*********************************************************/
  const double OneOverSqr2=0.7071067811865475;
  bool almout;
  double **gaus0, **gaus1;
  gsl_rng *rnd;
  int lmax, lmin;
  Alm<xcomplex <double> > *aflm;
  
  lmax = config.readi("LMAX");
  lmin = config.readi("LMIN");

  // Allocate memory:
  gaus0 = matrix<double>(0,Nfields-1, 0,1); // Complex random variables, [0] is real, [1] is imaginary part.
  gaus1 = matrix<double>(0,Nfields-1, 0,1); 
  aflm = vector<Alm<xcomplex <double> > >(0,Nfields-1); // Allocate Healpix Alm objects and set their size and initial value.
  for (i=0; i<Nfields; i++) {
    aflm[i].Set(lmax,lmax);
    for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) aflm[i](l,m).Set(0,0);
  }
  
  // - Set random number generator
  rnd = gsl_rng_alloc(gsl_rng_mt19937);
  if (rnd==NULL) error("corrlnfields: gsl_rng_alloc failed!");
  gsl_rng_set(rnd, config.readi("RNDSEED"));    // set random seed

  // LOOP over l's:
  j=0; // Will count number of Cholesky failures.
  for (l=lmin; l<=lmax; l++) {
    cout << "** Working with cov. matrix for l="<<l<<":\n";
    // Check if cov. matrix is positive definite and performs regularization if necessary:
    RegularizeCov(CovByl[l], config);
    // Output regularized matrix if requested:
    if (config.reads("REG_COVL_PREFIX")!="0") {
	filename=config.reads("REG_COVL_PREFIX")+"l"+ZeroPad(l,lmax)+".dat";
	GeneralOutput(CovByl[l], filename); 
    }

    // Perform a Cholesky decomposition:
    cout << "   Performing a Cholesky decomposition... "; cout.flush();
    status = gsl_linalg_cholesky_decomp(CovByl[l]);
    if (status==GSL_EDOM) { warning("Cholesky decomposition failed: matrix is not positive-definite."); j++; }
    cout << "done.\n";
    // Only proceed if Cholesky was successfull:
    if (status!=GSL_EDOM) { 
      // Output file if requested:
      if (config.reads("CHOLESKY_PREFIX")!="0") {
	filename=config.reads("CHOLESKY_PREFIX")+"l"+ZeroPad(l,lmax)+".dat";
	GeneralOutput(CovByl[l], filename); 
      }
      
      cout << "   Generating random auxiliary alm's...   "; cout.flush();
      // Generate m=0:
      m=0;
      for (i=0; i<Nfields; i++) {
	gaus0[i][0] = gsl_ran_gaussian(rnd, 1.0);
	gaus0[i][1] = 0.0;
      }
      CorrGauss(gaus1, CovByl[l], gaus0);
      for (i=0; i<Nfields; i++) aflm[i](l,m).Set(gaus1[i][0], gaus1[i][1]);
      // LOOP over m>0 for a fixed l:
      for (m=1; m<=l; m++) {
	// Generate independent 1sigma complex random variables:
	for (i=0; i<Nfields; i++) {
	  gaus0[i][0] = gsl_ran_gaussian(rnd, OneOverSqr2);
	  gaus0[i][1] = gsl_ran_gaussian(rnd, OneOverSqr2);
	}
	// Generate correlated complex gaussian variables according to CovMatrix:
	CorrGauss(gaus1, CovByl[l], gaus0);
	// Save alm to tensor:
	for (i=0; i<Nfields; i++) aflm[i](l,m).Set(gaus1[i][0], gaus1[i][1]);   
      } // End of LOOP over m's.
      cout << "done.\n";
    } // End of successfull Cholesky block. 
  } // End of LOOP over l's.
  free_matrix(gaus0,0,Nfields-1,0,1);
  free_matrix(gaus1,0,Nfields-1,0,1);
  free_GSLMatrixArray(CovByl, Nls);
  // Exit if any Cholesky failed:
  if (j>0) {sprintf(message,"Cholesky decomposition failed %d times.",j); error(message);} 

  // If requested, write alm's to file:
  GeneralOutput(aflm, config, "AUXALM_OUT", N1, N2);

  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="CHOLESKY_PREFIX" || 
      config.reads("EXIT_AT")=="AUXALM_OUT"      ||
      config.reads("EXIT_AT")=="REG_COVL_PREFIX" ) {
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
  cout << "Allocating memory for pixel maps...              "; cout.flush();
  nside   = config.readi("NSIDE");
  npixels = 12*nside*nside;
  mapf=vector<Healpix_Map<double> >(0,Nfields-1);
  for(i=0; i<Nfields; i++) mapf[i].SetNside(nside, RING); 		
  cout << "done.\n";
  // Generate maps from alm's for each field:
  cout << "Generating maps from alm's...                    "; cout.flush();
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
    cout << "LOGNORMAL realizations: exponentiating pixels... "; cout.flush();
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
  int Ngalaxies, gali, pixelNgal, PartialNgal, **catSet;
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
