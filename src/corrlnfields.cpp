/* corrlnfields: Written by Henrique S. Xavier on Nov-2014
   e-mail: hsxavier@if.usp.br
 */

#include <iostream>
#include "corrlnfields_aux.hpp" // Auxiliary functions made for this program.
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

/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl; using std::string; using std::ofstream; // Basic stuff.
  using namespace ParDef; ParameterList config;                             // Easy configuration file use.
  Cosmology cosmo;                                                          // Cosmological parameters.
  char message[100], message2[100];                                         // Handling warnings and errors.
  const int fgalaxies=1, fshear=2;                                          // Field type identification.
  std::string filename, tempstr;
  std::ofstream outfile;                                                    // File for output.
  std::ifstream infile;                                                     // File for input.
  enum simtype {gaussian, lognormal}; simtype dist;                         // For specifying simulation type.
  gsl_matrix *CovMatrix, **CovByl; 
  long CovSize;
  int status, i, j, l, m, Nfields, mmax, *ftype;
  double *means, *shifts, **aux, **zrange; 
  long long1, long2;
  FILE* stream; int NinputCls; std::string *filelist;                       // To list input Cls.
  gsl_set_error_handler_off();                                              // !!! All GSL return messages MUST be checked !!!


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
  // - Lognormal or Gausian realizations:
  if (config.reads("DIST")=="LOGNORMAL") dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN") dist=gaussian;
  else error("corrlnfields: unknown DIST: "+config.reads("DIST"));
  
  // Listing files to use based on CL_PREFIX:
  // Find out how many input C(l)'s we have.
  sprintf(message, "ls %s* | wc -l", config.reads("CL_PREFIX").c_str()); 
  stream = popen(message, "r");
  if ((stream=popen(message, "r")) == NULL) error("corrlnfields: cannot 'ls | wc' output");
  fscanf(stream, "%d", &NinputCls); pclose(stream);
  // Allocate memory:
  filelist = vector<std::string>(0,NinputCls-1);
  // Get list of input C(l) files:
  sprintf(message, "ls %s*", config.reads("CL_PREFIX").c_str());
  if ((stream=popen(message, "r")) == NULL) error("corrlnfields: cannot pipe 'ls' output");
  for (i=0; i<NinputCls; i++) {
    fscanf(stream, "%s", message);
    filelist[i].assign(message);
  }
  pclose(stream);

  
  /********************************************/
  /*** PART 1: Load C(l)s and organize them ***/
  /********************************************/
  void CountEntries(std::string filename, long *nr, long *nc);
  void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2);
  int a1, a2, b1, b2, N1, N2, Nlinput, **fnz, **NentMat;
  long *Nentries, ncols;
  double ***ll, ***Cov, *wrapper[2];
  bool *fnzSet, **IsSet;
  
  // Get file list and find out how many C(l)s there are:  
  N1=0; N2=0, Nlinput=0;
  Nentries = vector<long>(0,NinputCls-1);
  for (i=0; i<NinputCls; i++) {
    getcovid(filelist[i], &a1, &a2, &b1, &b2);
    if (a1>N1) N1=a1; if (b1>N1) N1=b1;                // Get number of fields.
    if (a2>N2) N2=a2; if (b2>N2) N2=b2;                // Get number of z bins.
    CountEntries(filelist[i], &(Nentries[i]), &ncols); // Get number of Nls.
    if (ncols!=2) error("corrlnfields: wrong number of columns in file "+filename);
    if (Nentries[i]>Nlinput) Nlinput=Nentries[i];          // Record maximum number of ls.
  }
  cout << "Nfields: " << N1 << " Nzs: " << N2 << endl;
  
  // Allocate memory to store C(l)s:
  // First two indexes are CovMatrix indexes and last is for ll.
  // fnz stores the order that the fields are stored in CovMatrix.
  fnz     =     matrix<int>(0, N1*N2-1, 0, 1);                 // Records what field is stored in each element of CovMatrix.
  fnzSet  =    vector<bool>(0, N1*N2-1);                       // For bookkeeping.
  ll      = tensor3<double>(0, N1*N2-1, 0, N1*N2-1, 0, Nlinput); // Records the ll for each C(l) file. 
  Cov     = tensor3<double>(0, N1*N2-1, 0, N1*N2-1, 0, Nlinput); // Records the C(l) for each C(l) file.
  IsSet   =    matrix<bool>(0, N1*N2-1, 0, N1*N2-1);           // For bookkeeping.
  NentMat =     matrix<int>(0, N1*N2-1, 0, N1*N2-1);           // Number of C(l) entries in file.
  for(i=0; i<N1*N2; i++) for(j=0; j<N1*N2; j++) IsSet[i][j]=0;
  for(i=0; i<N1*N2; i++) fnzSet[i]=0;
  
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
  for (i=0; i<N1*N2; i++) if (fnzSet[i]==0) error("corrlnfields: some position in CovMatrix is unclaimed.");
  free_vector(fnzSet, 0, N1*N2-1);
  // If positions are OK and output required, print them out:
  if (config.reads("FLIST_OUT")!="0") {
    outfile.open(config.reads("FLIST_OUT").c_str());
    if (!outfile.is_open()) error("corrlnfields: cannot open FLIST_OUT file.");
    PrintTable(fnz, N1*N2, 2, &outfile);
    outfile.close();
    cout << "Written field list to "+config.reads("FLIST_OUT")<<endl;
  }
  

  /**************************************************/
  /*** PART 2: Prepare for Cholesky decomposition ***/
  /**************************************************/
  double *tempCl, *LegendreP, *workspace, *xi, lsup, supindex, *theta, *DLTweights, *lls;
  const int HWMAXL = 10000000; int lastl = HWMAXL, Nls;
  
  // Load means, shifts, type and z range data file:
  cout << "Loading means and shifts from file "+config.reads("FIELDS_INFO")+":\n";
  aux   = LoadTable<double>(config.reads("FIELDS_INFO"), &long1, &long2); 
  Nfields = (int)long1; // From now on will use Nfields instead of N1*N2, so the check below is important!
  if (Nfields != N1*N2) error("corrlnfields: number of means and shifts do not match number of C(l)s.");
  fnzSet = vector<bool>(0, Nfields-1); for (i=0; i<Nfields; i++) fnzSet[i]=0;
  means  = vector<double>(0, Nfields-1);
  ftype  = vector<int>(0, Nfields-1);
  zrange = matrix<double>(0,Nfields-1, 0,1);
  if (dist==lognormal) shifts = vector<double>(0, Nfields-1);
  for (j=0; j<Nfields; j++) {
    fz2n((int)aux[j][0], (int)aux[j][1], &i, N1, N2); // Find conventional position of field in arrays.
    if (fnzSet[i]==1) error ("corrlnfields: found more than one mean & shift entry for the same f-z.");
    fnzSet[i] = 1; 
    means[i]  = aux[j][2];  
    if (dist==lognormal) shifts[i] = aux[j][3];
    ftype[i] = (int)aux[j][4];
    zrange[i][0] = aux[j][5]; zrange[i][1] = aux[j][6]; 
  }
  for (i=0; i<Nfields; i++) if (fnzSet[i]!=1) error("corrlnfields: the properties of a field were not set.");
  for (i=0; i<Nfields; i++) if (zrange[i][0]>zrange[i][1]) error("corrlnfields: zmin > zmax for a field.");
  if (dist==lognormal) for (i=0; i<Nfields; i++) if(means[i]+shifts[i]<=0) { // Sanity check.
	printf(message, "corrlnfields: mean+shift at position %d must be greater than zero.", i); error(message);
      }
  free_vector(fnzSet, 0, Nfields-1);
  free_matrix(aux, 0, Nfields-1, 0, long2-1);
  cout << "Done.\n";
  
  // Look for the maximum l value described by all C(l)s:
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]==1) {
	if (ll[i][j][NentMat[i][j]-1]>HWMAXL) error ("corrlnfields: too high l in C(l)s: increase HWMAXL.");
	if (ll[i][j][NentMat[i][j]-1]<lastl) lastl = (int)ll[i][j][NentMat[i][j]-1];
      }
  Nls=lastl+1; // l=0 is needed for DLT. Nls is known as 'bandwidth' (bw) in s2kit 1.0 code.

  // Allocate gsl_matrices that will receive covariance matrices for each l.
  cout << "Allocating data-cube necessary for Cholesky decomposition... "; cout.flush();
  tempCl = vector<double>(0, lastl);
  CovByl = GSLMatrixArray(Nls, Nfields, Nfields);
  cout << "done.\n";


  /*****************************************************************/
  /*** PART 3: Compute auxiliary gaussian C(l)s if LOGNORMAL     ***/
  /*****************************************************************/
  if (dist==lognormal) {
    cout << "LOGNORMAL realizations: will compute auxiliary gaussian C(l)s:\n";
    // Loads necessary memory:
    cout << "Allocating extra memory... "; cout.flush();
    workspace  = vector<double>(0, 16*Nls-1);
    LegendreP  = vector<double>(0, 2*Nls*Nls-1);
    xi         = vector<double>(0, 2*Nls-1);
    theta      = vector<double>(0, 2*Nls-1);
    lls        = vector<double>(0, lastl);
    DLTweights = vector<double>(0, 4*Nls-1);
    // Initialize vectors:
    for (i=0; i<=lastl; i++) lls[i]=(double)i;
    ArcCosEvalPts(2*Nls, theta);
    for (i=0; i<2*Nls; i++) theta[i] = theta[i]*180.0/M_PI; 
    cout << "done.\n";
    // Loads C(l) exponential suppression:
    lsup     = config.readd("SUPPRESS_L");
    supindex = config.readd("SUP_INDEX"); 
    // Load s2kit 1.0 Legendre Polynomials:
    cout << "Generating table of Legendre polynomials and sampling angles... "; cout.flush();
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
	GetAllLs(ll[i][j], Cov[i][j], NentMat[i][j], tempCl, lastl);
	cout << "              done.\n";
	
	if (dist==lognormal) {              /** LOGNORMAL ONLY **/
	  // Compute correlation function Xi(theta):
	  cout << "   DLT (inverse) to obtain the correlation function... "; cout.flush();
	  ModCl4DLT(tempCl, lastl, lsup, supindex);
	  Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	  cout << "  done.\n";
	  if (config.reads("XIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("XIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    cout << "   Correlation function written to "+filename<<endl;
	  }
	  // Transform Xi(theta) to auxiliary gaussian Xi(theta):
	  cout << "   Computing associated gaussian correlation function... "; cout.flush(); 
	  status=GetGaussCorr(xi, xi, 2*Nls, means[i], shifts[i], means[j], shifts[j]);
	  cout << "done.\n";
	  if (status==EDOM) error("corrlnfields: GetGaussCorr found bad log arguments.");
	  if (config.reads("GXIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GXIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    cout << "   Associated Gaussian correlation function written to "+filename<<endl;
	  }
	  // Transform Xi(theta) back to C(l):
	  cout << "   DLT (forward) to obtain the angular power spectrum... "; cout.flush(); 
	  Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	  ApplyClFactors(tempCl, Nls);
	  cout << "done.\n";
	  if (config.reads("GCLOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GCLOUT_PREFIX"), i, j, N1, N2, lls, tempCl, Nls);
	    cout << "   C(l) for auxiliary Gaussian variables written to "+filename<<endl;
	  }	  
	}                                 /** END OF LOGNORMAL ONLY **/ 
	
	// Save auxiliary C(l):
	for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j]=tempCl[l];		
      } // End of LOOP over C(l)[i,j] that were set.
  
  // Freeing memory: from now on we only need CovByl, means, shifts, fnz.
  cout << "Massive memory deallocation... "; cout.flush();
  free_vector(tempCl, 0, lastl);
  free_tensor3(Cov,    0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_tensor3(ll,     0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_matrix(NentMat, 0, Nfields-1, 0, Nfields-1);
  if (dist==lognormal) {
    free_vector(workspace, 0, 16*Nls-1);
    free_vector(LegendreP, 0, 2*Nls*Nls-1);
    free_vector(xi, 0, 2*Nls-1);
    free_vector(theta, 0, 2*Nls-1);
    free_vector(lls, 0, lastl);
    free_vector(DLTweights, 0, 4*Nls-1); 
  }
  cout << "done.\n";

  
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
  for (l=lmin; l<=lmax; l++) {
    cout << "** Working with cov. matrix for l="<<l<<":\n";
    
    // Perform a Cholesky decomposition:
    cout << "   Performing a Cholesky decomposition... "; cout.flush();
    status = gsl_linalg_cholesky_decomp(CovByl[l]);
    if (status==GSL_EDOM) error("Cholesky decomposition failed: matrix is not positive-definite.");
    cout << "done.\n";
    // Output file if requested:
    if (config.reads("CHOLESKY_PREFIX")!="0") {
      filename=config.reads("CHOLESKY_PREFIX")+"l"+ZeroPad(l,lmax)+".dat";
      outfile.open(filename.c_str());
      if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
      else { 
	PrintGSLMatrix(CovByl[l], &outfile); outfile.close();
	cout << "   Cholesky decomposition output matrix written to " << filename << endl;
      }  
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
    
  } // End of LOOP over l's.
  free_matrix(gaus0,0,CovSize-1,0,1);
  free_matrix(gaus1,0,CovSize-1,0,1);
  free_GSLMatrixArray(CovByl, Nls);

  // If requested, write alm's to file:
  if (config.reads("AUXALM_OUT")!="0") {
    filename = config.reads("AUXALM_OUT");
    outfile.open(filename.c_str());
    if (!outfile.is_open()) 
      warning("corrlnfields: cannot open "+filename+" file.");
    outfile << SampleHeader(fnz, Nfields) <<endl<<endl;
    lmin = config.readi("LRANGE_OUT", 0);
    lmax = config.readi("LRANGE_OUT", 1);
    mmax = config.readi("MMAX_OUT");
    if (mmax>lmin) error ("corrlnfields: current code only allows MMAX_OUT <= LMIN_OUT.");
    // Output all alm's:
    if (mmax<0) {
      for(l=lmin; l<=lmax; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< aflm[i](l,m).re<<" "<<std::setprecision(10)<< aflm[i](l,m).im;
	  outfile<<endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lmin; l<=lmax; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< aflm[i](l,m).re<<" "<<std::setprecision(10)<< aflm[i](l,m).im;
	  outfile<<endl;
	}  
    }
    outfile.close();
    cout << "Auxiliary alm's written to "+filename<<endl;
  }
  
  /******************************/
  /*** Part 5: Map generation ***/
  /******************************/
  int nside, npixels;
  Healpix_Map<double> *mapf;
  double expmu, gmean, gvar;
  pointing coord;
  int field, z;
  char *arg[5];
  char opt1[]="-bar", val1[]="1";

  // Allocate memory for pixel maps:
  cout << "Allocating memory for pixel maps... "; cout.flush();
  nside   = config.readi("NSIDE");
  npixels = 12*nside*nside;
  mapf=vector<Healpix_Map<double> >(0,Nfields-1);
  for(i=0; i<Nfields; i++) mapf[i].SetNside(nside, RING); 		
  cout << "done.\n";
  // Generate maps from alm's for each field:
  cout << "Generating maps from alm's... "; cout.flush();
  for(i=0; i<Nfields; i++) alm2map(aflm[i],mapf[i]);
  cout << "done.\n";

  // Write auxiliary map to file as a table if requested:
  if (config.reads("AUXMAP_OUT")!="0") {
    filename = config.reads("AUXMAP_OUT");
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
    else {
      outfile << "# theta, phi";
      for (i=0; i<Nfields; i++) {
	n2fz(i, &field, &z, N1, N2);
	outfile << ", f"<<field<<"z"<<z;
      }
      outfile << endl;
      for (j=0; j<npixels; j++) {
	coord = mapf[0].pix2ang(j);
	outfile << coord.theta <<" "<< coord.phi;
	for (i=0; i<Nfields; i++) outfile <<" "<< mapf[i][j];
	outfile << endl;
      }
      outfile.close();
      cout << "Auxiliary catalog written to " << filename << endl;
    }  
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
  if (config.reads("MAP_OUT")!="0") {
    filename = config.reads("MAP_OUT");
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
    else {
      outfile << "# theta, phi";
      for (i=0; i<Nfields; i++) {
	n2fz(i, &field, &z, N1, N2);
	outfile << ", f"<<field<<"z"<<z;
      }
      outfile << endl;
      for (j=0; j<npixels; j++) {
	coord = mapf[0].pix2ang(j);
	outfile << coord.theta <<" "<< coord.phi;
	for (i=0; i<Nfields; i++) outfile <<" "<< mapf[i][j];
	outfile << endl;
      }
      outfile.close();
      cout << "Catalog written to " << filename << endl;
    }  
  }
  
  // Map output to fits and/or tga files:
  if (config.reads("MAPFITS_PREFIX")!="0") {
    // Write to FITS:
    cout << "Writing maps to fits files:\n";
    tempstr  = config.reads("MAPFITS_PREFIX");
    for (i=0; i<Nfields; i++) {
      sprintf(message, "%sf%dz%d.fits", tempstr.c_str(), fnz[i][0], fnz[i][1]);
      filename.assign(message);
      sprintf(message2, "rm -f %s", message);
      system(message2); // Have to delete previous fits files first.
      write_Healpix_map_to_fits(filename,mapf[i],planckType<double>());
      cout << "Map for field ["<<i<<"] written to "<<filename<<endl;
      // Write to TGA if requested:
      if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
	sprintf(message2, "%sf%dz%d.tga", tempstr.c_str(), fnz[i][0], fnz[i][1]);
	arg[1]=message; arg[2]=message2; arg[3]=opt1; arg[4]=val1;
	map2tga_module(4, (const char **)arg);
	cout << "Map for field ["<<i<<"] written to "<<message2<<endl;
	if (config.readi("FITS2TGA")==2) {
	  sprintf(message2, "rm -f %s", message);
	  system(message2);
	  cout << "Deleted "<<filename<<endl;;
	}
      }
    }
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
    filename = config.reads("RECOVALM_OUT");  
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("corrlnfields: cannot open "+filename+" file.");
    outfile << SampleHeader(fnz, Nfields) <<endl<<endl;
    lmin = config.readi("LRANGE_OUT", 0);
    lmax = config.readi("LRANGE_OUT", 1);
    mmax = config.readi("MMAX_OUT");
    if (mmax>lmin) error ("corrlnfields: current code only allows MMAX_OUT <= LMIN_OUT.");
    // Output all alm's:
    if (mmax<0) {
      for(l=lmin; l<=lmax; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< aflm[i](l,m).re<<" "<<std::setprecision(10)<< aflm[i](l,m).im;
	  outfile<<endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lmin; l<=lmax; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< aflm[i](l,m).re<<" "<<std::setprecision(10)<< aflm[i](l,m).im;
	  outfile<<endl;
	}  
    } 
    outfile.close();
    cout << "Recovered alm's written to "+filename<<endl;
  }
  free_vector(aflm, 0, Nfields-1);


  /**********************************/
  /*** Part 6: catalog generation ***/
  /**********************************/
  Healpix_Map<double> *selection;
  double PixelSolidAngle=12.56637061435917/npixels; // 4pi/npixels.

  // Read selection functions from FITS files:
  cout << "Reading selection functions from files... "; cout.flush();
  selection = vector<Healpix_Map<double> >(0,Nfields-1);
  tempstr   = config.reads("SELEC_PREFIX");
  for (i=0; i<Nfields; i++) {
    sprintf(message, "%sf%dz%d.fits", tempstr.c_str(), fnz[i][0], fnz[i][1]);
    filename.assign(message);
    read_Healpix_map_from_fits(filename, selection[i]);
    if (selection[i].Nside()!=mapf[i].Nside()) 
      error("corrlnfields: "+filename+" does not have the same Nside as the corresponding map.");
    if (selection[i].Scheme()!=mapf[i].Scheme())
      error("corrlnfields: "+filename+" does not have the same ordering scheme as the corresponding map.");
  }
  cout << "done.\n";
  // If SELEC_TYPE==FRACTION, multiply it by the mean projected density:
  if (config.reads("SELEC_TYPE")=="FRACTION") error ("corrlnfields: SELEC_TYPE FRACTION not implemented yet.");
  else if (config.reads("SELEC_TYPE")!="DENSITY") error ("corrlnfields: unknown SELEC_TYPE option.");
  
  // Poisson Sampling the galaxy fields:
  if (config.readi("POISSON")==1) {
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	cout << "Poisson sampling f"<<fnz[i][0]<<"z"<<fnz[i][1]<<"... "; cout.flush();
	for(j=0; j<npixels; j++) mapf[i][j] = gsl_ran_poisson(rnd, selection[i][j]*(1.0+mapf[i][j])/**PixelSolidAngle*/);
	cout << "done.\n";
      }
  }
  // Just output the expected number density, if requested:
  else if (config.readi("POISSON")==0) {
    for (i=0; i<Nfields; i++) if (ftype[i]==fgalaxies) {
	cout << "Using expected number density for f"<<fnz[i][0]<<"z"<<fnz[i][1]<<"... "; cout.flush();
	for(j=0; j<npixels; j++) mapf[i][j] = selection[i][j]*(1.0+mapf[i][j])/**PixelSolidAngle*/;
	cout << "done.\n";
      }
  }
  else error ("corrlnfields: unknown POISSON option.");
  
  cout << "Testing pixel boundaries:\n";
  //std::vector<vec3> corner;
  double Theta, Phi;
  pointing rndang;
  Healpix_Map<double> BaseMap, FineMap;
  int Nsamples=1000;
  BaseMap.SetNside(2, RING);
  BaseMap.fill(0);
  FineMap.SetNside(64, RING);
  FineMap.fill(0);
  i=47; // Escolhemos um pixel do BaseMap e o destacamos.
  BaseMap[i]=1.0;
  for (j=0; j<Nsamples*64*64; j++) {
    rndang = RandAngInPix(rnd, BaseMap, i);
    l = FineMap.ang2pix(rndang);
    FineMap[l] = FineMap[l] + 1.0/Nsamples;
  }
  write_Healpix_map_to_fits("base.fits",BaseMap,planckType<double>());
  write_Healpix_map_to_fits("fine.fits",FineMap,planckType<double>());
  
  

  // End of the program
  free_matrix(fnz, 0,N1*N2-1, 0,1);
  free_vector(ftype, 0, Nfields-1);
  free_matrix(zrange,0, Nfields-1,0,1);
  free_vector(mapf, 0, Nfields-1);
  gsl_rng_free(rnd);
  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}
