#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling and memory allocation.
#include "s2kit10_naive.hpp"    // For Discrete Legendre Transforms.
#include <gsl/gsl_matrix.h>     // gsl_matrix.
#include <math.h>               // exp, log...
#include "GeneralOutput.hpp"
#include "RegularizeCov.hpp"
#include "corrlnfields_aux.hpp" // For n2fz functions, Maximum, etc..


/*** Transforms a correlation function of gaussian variables gXi into a corr. function of corresponding lognormal variables lnXi ***/
void GetLNCorr(double *lnXi, double *gXi, int XiLength, double mean1, double shift1, double mean2, double shift2) {
  int i;
  
  for (i=0; i<XiLength; i++) lnXi[i] = ( exp(gXi[i]) - 1.0 ) * (mean1+shift1) * (mean2+shift2);
}


/*** Transforms a correlation function of lognormal variables lnXi into a corr. function of associated gaussian variables gXi ***/
int GetGaussCorr(double *gXi, double *lnXi, int XiLength, double mean1, double shift1, double mean2, double shift2) {
  int i, status=0;
  double arg, bad=-666.0;
  char message[100];
  
  for (i=0; i<XiLength; i++) {
    arg = 1.0 + lnXi[i]/(mean1+shift1)/(mean2+shift2);
    if (arg <= 0) {
      sprintf(message, "GetGaussCorr: lnXi[%d] leads to bad log argument, gXi[%d] set to %g.", i, i, bad);
      warning(message);
      status=EDOM;
      gXi[i] = bad;
    }
    else gXi[i] = log(arg);
  }
  return status;
}


/*** Get four numbers separated by characters that specify the fields and redshifts
     of the correlation function. ***/
void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2) {
  int i=0, num, index, fileL;
  fileL=filename.length();
  // LOOP over the four indexes that indentifies the C(l):
  for(index=1; index<=4; index++) {
    num=0;
    // Find a number:
    while (isdigit(filename.c_str()[i])==0) {i++; if(i>=fileL) error("getcovid: cannot find four numbers.");}
    // Read the number:
    while (isdigit(filename.c_str()[i])!=0) {num = num*10 + (filename.c_str()[i]-'0'); i++;}
    // Save it as an index:
    switch (index) {
    case 1: *a1 = num; break;
    case 2: *a2 = num; break;
    case 3: *b1 = num; break;
    case 4: *b2 = num; break;
    }
  }
  // Check if there are more numbers in filename:
  while (i<=fileL) {
    if (isdigit(filename.c_str()[i])!=0) error("getcovid: found more numbers than expected.");
    i++;
  }
}


/*** Find out number of columns and rows in file ***/
void CountEntries(std::string filename, long *nr, long *nc) {
  using std::ifstream;
  using std::string;
  using std::istringstream;
  using std::ostringstream;
  long nrows=0, ncols=0;
  ifstream file;
  istringstream inputline; ostringstream outputline;
  string word, phrase;
  
  // Open file
  file.open(filename.c_str());
  if (!file.is_open()) error("CountEntries: cannot open file.");
  
  // Count lines and columns:
  getline(file,phrase);
  outputline << phrase;
  inputline.str(outputline.str());
  while (inputline >> word) ncols++;
  while(!file.eof()) {getline(file,phrase); if (phrase.length()>0) nrows++;}

  file.close();
  *nr=nrows+1;
  *nc=ncols;
}


/*** Export function y(x) for the field combination [i,j] to file ***/
std::string PrintOut(std::string prefix, int i, int j, int N1, int N2, double *x, double *y, int length) {
  int a1, a2, b1, b2;
  char message[100];
  std::string filename;
  double *wrapper[2];
  std::ofstream outfile;

  wrapper[0] =  x;
  wrapper[1] =  y;

  n2fz(i, &a1, &a2, N1, N2); n2fz(j, &b1, &b2, N1, N2);
  sprintf(message, "%sf%dz%df%dz%d.dat", prefix.c_str(),a1,a2,b1,b2);
  filename.assign(message);

  outfile.open(message);
  if (!outfile.is_open()) error("PrintOut: cannot open file "+filename);
  PrintVecs(wrapper, length, 2, &outfile);
  outfile.close();

  return filename;
}


/*** Wraps all processing of the input Cls up to the gaussian covariance matrices for each l ***/
int ClProcess(gsl_matrix ***CovBylAddr, double *means, double *shifts, int N1, int N2, int *NlsOut, const ParameterList & config) {
  using namespace definitions;                          // Global definitions.
  using std::cout; using std::endl;                     // Basic stuff.
  simtype dist;                                         // For specifying simulation type.
  char message[200];                                    // Writing outputs with sprintf.
  std::ofstream outfile;                                // File for output.
  FILE* stream; int NinputCls; std::string *filelist;   // To list input Cls.
  int i, j, l, m, status, Nfields, Nls;
  std::string filename;
  bool *fnzSet;
  gsl_matrix **CovByl;

  // Getting general information:
  Nfields=N1*N2;
  if (config.reads("DIST")=="LOGNORMAL") dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN") dist=gaussian;


  /********************************************/
  /*** PART 1: Load C(l)s and organize them ***/
  /********************************************/
  const int HWMAXL = 10000000; int lastl = HWMAXL;
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
  if (config.reads("EXIT_AT")=="FLIST_OUT") return 1;
 

  // Look for the maximum l value described by all C(l)s:
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]==1) {
	if (ll[i][j][NentMat[i][j]-1]>HWMAXL) error ("corrlnfields: too high l in C(l)s: increase HWMAXL.");
	if (ll[i][j][NentMat[i][j]-1]<lastl) lastl = (int)ll[i][j][NentMat[i][j]-1];
      }
  Nls=lastl+1; // l=0 is needed for DLT. Nls is known as 'bandwidth' (bw) in s2kit 1.0 code.
  (*NlsOut)=Nls;

  
  // Allocate gsl_matrices that will receive covariance matrices for each l.
  cout << "Allocating data-cube needed for Cholesky decomposition... "; cout.flush();
  CovByl = GSLMatrixArray(Nls, Nfields, Nfields);
  (*CovBylAddr) = CovByl;
  cout << "done.\n";
 
  /*****************************************************************/
  /*** PART 2: Compute auxiliary gaussian C(l)s if LOGNORMAL     ***/
  /*****************************************************************/
  double *tempCl, *LegendreP, *workspace, *xi, lsup, supindex, *theta, *DLTweights, *lls;

  if (dist==lognormal) {
    cout << "LOGNORMAL realizations: will compute auxiliary gaussian C(l)s:\n";
    // Loads necessary memory:
    cout << "Allocating memory for DLT...                              "; cout.flush();
    lls        = vector<double>(0, lastl);
    workspace  = vector<double>(0, 16*Nls-1);
    LegendreP  = vector<double>(0, 2*Nls*Nls-1);
    DLTweights = vector<double>(0, 4*Nls-1);
    // Initialize vectors:
    for (i=0; i<=lastl; i++) lls[i]=(double)i;
    cout << "done.\n";
    // angle theta is only necessary for output:
    if (config.reads("XIOUT_PREFIX")!="0" || config.reads("GXIOUT_PREFIX")!="0") {
      cout << "Generating table of sampling angles...                    "; cout.flush();
      theta    = vector<double>(0, 2*Nls-1);
      ArcCosEvalPts(2*Nls, theta);
      for (i=0; i<2*Nls; i++) theta[i] = theta[i]*180.0/M_PI;
      cout << "done.\n";
    } 
    
    // Loads C(l) exponential suppression:
    lsup     = config.readd("SUPPRESS_L");
    supindex = config.readd("SUP_INDEX"); 
    // Load s2kit 1.0 Legendre Polynomials:
    cout << "Generating table of Legendre polynomials...               "; cout.flush();
    PmlTableGen(Nls, 0, LegendreP, workspace);
    free_vector(workspace, 0, 16*Nls-1);
    cout << "done.\n";
    // Compute s2kit 1.0 Discrete Legendre Transform weights:
    cout << "Calculating forward DLT weights...                        "; cout.flush();
    makeweights(Nls, DLTweights);
    cout << "done.\n";
  }

  // LOOP over all C(l)s already set.
  if (dist==lognormal) {cout << "Transforming C(l)s for the auxiliary Gaussian ones...     "; cout.flush();}
  else cout << "Interpolating C(l)s for all l's... "; cout.flush();
#pragma omp parallel for collapse(2) schedule(dynamic) private(tempCl, xi, workspace, filename, l)
  for(i=0; i<Nfields; i++)
    for(j=0; j<Nfields; j++) 
      if (IsSet[i][j]==1) {
 
	//cout << "** Transforming C(l) in ["<<i<<", "<<j<<"]:\n";
	
       	// Temporary memory allocation:
	tempCl    = vector<double>(0, lastl);
	xi        = vector<double>(0, 2*Nls-1);
	workspace = vector<double>(0, 2*Nls-1);

	// Interpolate C(l) for every l; input C(l) might not be like that:
	//cout << "   Interpolating input C(l) for all l's...  "; cout.flush();
	GetAllLs(ll[i][j], Cov[i][j], NentMat[i][j], tempCl, lastl, config.readi("EXTRAP_DIPOLE"));
	//cout << "              done.\n";
	
	if (dist==lognormal) {              /** LOGNORMAL ONLY **/

	  // Compute correlation function Xi(theta):
	  //cout << "   DLT (inverse) to obtain the correlation function...  "; cout.flush();
	  ModCl4DLT(tempCl, lastl, lsup, supindex);
	  Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	  //cout << "  done.\n";
	  if (config.reads("XIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("XIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    //cout << ">> Correlation function written to "+filename<<endl;
	  }

	  // Transform Xi(theta) to auxiliary gaussian Xi(theta):
	  //cout << "   Computing associated gaussian correlation function...  "; cout.flush(); 
	  status=GetGaussCorr(xi, xi, 2*Nls, means[i], shifts[i], means[j], shifts[j]);
	  //cout << "done.\n";
	  if (status==EDOM) error("corrlnfields: GetGaussCorr found bad log arguments.");
	  if (i==j && xi[0]<0) warning("corrlnfields: auxiliary field variance is negative.");
	  if (config.reads("GXIOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GXIOUT_PREFIX"), i, j, N1, N2, theta, xi, 2*Nls);
	    //cout << ">> Associated Gaussian correlation function written to "+filename<<endl;
	  }

	  // Transform Xi(theta) back to C(l):
	  //cout << "   DLT (forward) to obtain the angular power spectrum...  "; cout.flush(); 
	  Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	  ApplyClFactors(tempCl, Nls);
	  //cout << "done.\n";
	  if (config.reads("GCLOUT_PREFIX")!="0") { // Write it out if requested:
	    filename=PrintOut(config.reads("GCLOUT_PREFIX"), i, j, N1, N2, lls, tempCl, Nls);
	    //cout << ">> C(l) for auxiliary Gaussian variables written to "+filename<<endl;
	  }	  
	}                                 /** END OF LOGNORMAL ONLY **/ 
	
	// Save auxiliary C(l):
	for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j]=tempCl[l];
	
	// Temporary memory deallocation:
	free_vector(tempCl, 0, lastl);
	free_vector(xi, 0, 2*Nls-1);
	free_vector(workspace, 0, 2*Nls-1);
      } // End of LOOP over C(l)[i,j] that were set.
  cout << "done.\n";
  
  // Memory deallocation:
  free_tensor3(Cov,    0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_tensor3(ll,     0, Nfields-1, 0, Nfields-1, 0, Nlinput); 
  free_matrix(NentMat, 0, Nfields-1, 0, Nfields-1);
  if (config.reads("XIOUT_PREFIX")!="0" || config.reads("GXIOUT_PREFIX")!="0") free_vector(theta, 0, 2*Nls-1);

  // Output information:
  if (config.reads("XIOUT_PREFIX")!="0") 
    cout << ">> Correlation functions written to prefix "+config.reads("XIOUT_PREFIX")<<endl;
  if (config.reads("GXIOUT_PREFIX")!="0") 
    cout << ">> Associated Gaussian correlation functions written to prefix "+config.reads("GXIOUT_PREFIX")<<endl;
  if (config.reads("GCLOUT_PREFIX")!="0") 
    cout << ">> C(l)s for auxiliary Gaussian variables written to prefix "+config.reads("GCLOUT_PREFIX")<<endl;
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="XIOUT_PREFIX"  || 
      config.reads("EXIT_AT")=="GXIOUT_PREFIX" || 
      config.reads("EXIT_AT")=="GCLOUT_PREFIX") return 1;
    
  // Set Cov(l)[i,j] = Cov(l)[j,i]
  cout << "Set remaining cov. matrices elements based on symmetry... "; cout.flush(); 
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
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]!=1) {
	sprintf(message,"corrlnfields: [%d,%d] was not set.",i,j);
	error(message);
      }
  cout << "done.\n";
  free_matrix(IsSet, 0, Nfields-1, 0, Nfields-1);
  
  // Output covariance matrices for each l if requested:
  GeneralOutput(CovByl, config, "COVL_PREFIX", 0);
  if (config.reads("COVL_PREFIX")!="0") 
    cout << ">> Cov. matrices written to prefix "+config.reads("COVL_PREFIX")<<endl;
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="COVL_PREFIX") return 1;


  /****************************************************************************/
  /*** PART 3: Regularize (make them positive definite) covariance matrices ***/
  /****************************************************************************/
  gsl_matrix *gslm;
  double *MaxChange, MMax;
  int lmax, lmin, lMMax, lstart, lend;

  lmax = config.readi("LMAX");
  lmin = config.readi("LMIN");

  // If producing the regularized lognormal C(l)s, all ls must be regularized (skip l=0 because is set to zero).
  // Else, only the cov. matrices for requested ls are regularized.
  // Note that a very high exponential suppression of C(l)s make difficult to regularize matrices.
  if (dist==lognormal && config.reads("REG_CL_PREFIX")!="0") { lstart = 1;  lend = Nls-1; }
  else { lstart = lmin;  lend = lmax; }
  MaxChange = vector<double>(lstart, lend);
  
  cout << "Regularizing cov. matrices...                             "; cout.flush();

#pragma omp parallel for schedule(dynamic) private(gslm, filename)  
  for (l=lstart; l<=lend; l++) {

    // Check pos. defness, regularize if necessary, keep track of changes:
    gslm = gsl_matrix_alloc(Nfields, Nfields);
    gsl_matrix_memcpy(gslm, CovByl[l]);
    status = RegularizeCov(CovByl[l], config);
    MaxChange[l] = MaxFracDiff(CovByl[l], gslm);
    if (status==9) { 
      sprintf(message, "ClProcess: RegularizeCov for l=%d reached REG_MAXSTEPS with Max. change of %g.",l,MaxChange[l]); 
      warning(message);
    }
    gsl_matrix_free(gslm);
    // Output regularized matrix if requested:
    if (config.reads("REG_COVL_PREFIX")!="0") {
      filename=config.reads("REG_COVL_PREFIX")+"l"+ZeroPad(l,lend)+".dat";
      GeneralOutput(CovByl[l], filename, 0);
    }
  }
  cout << "done.\n";
 
  // Dump changes in cov. matrices to the screen:
  for (l=lmin; l<=lmax; l++) {
    MMax  = 0.0;
    lMMax = 0;
    if (MaxChange[l]>MMax) {MMax = MaxChange[l]; lMMax = l;}
  }
  cout << "Max. % change for "<<lmin<<"<=l<="<<lmax<<" at l="<<lMMax<<": "<<MMax<<endl;  
  free_vector(MaxChange, lstart, lend);
  // Output regularized matrices if requested:
  if (config.reads("REG_COVL_PREFIX")!="0") 
    cout << ">> Regularized cov. matrices written to prefix "+config.reads("REG_COVL_PREFIX")<<endl;
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="REG_COVL_PREFIX") return 1;

  /***********************************************************/
  /*** PART 4: Obtain regularized input Cls if requested   ***/
  /***********************************************************/

  if (config.reads("REG_CL_PREFIX")!="0") {
    if(dist==lognormal) {
      cout << "Computing regularized lognormal Cls...                    "; cout.flush();
      
      
      // LOOP over fields:
#pragma omp parallel for collapse(2) schedule(dynamic) private(tempCl, xi, workspace, filename, l)
      for(i=0; i<Nfields; i++)
	for(j=i; j<Nfields; j++) {
	  //cout << "** Transforming C(l) in ["<<i<<", "<<j<<"]:\n";

	  // Temporary memory allocation:
	  tempCl    = vector<double>(0, lastl);
	  xi        = vector<double>(0, 2*Nls-1);
	  workspace = vector<double>(0, 2*Nls-1);

	  // Copy the Cl to a vector:
	  for (l=0; l<Nls; l++) tempCl[l] = CovByl[l]->data[i*Nfields+j]; // tudo certo.
	  // Compute correlation function Xi(theta):
	  //cout << "   DLT (inverse) to obtain the correlation function...     "; cout.flush();
	  ModCl4DLT(tempCl, lastl, -1, -1); // Suppression not needed (it was already suppressed).
	  Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	  //cout << "done.\n";
	  // Get Xi(theta) for lognormal variables:
	  //cout << "   Getting correlation function for lognormal variables... "; cout.flush();
	  GetLNCorr(xi, xi, 2*Nls, means[i], shifts[i], means[j], shifts[j]);
	  //cout << "done.\n";
	  // Compute the Cls:
	  //cout << "   DLT (forward) to obtain the angular power spectrum...   "; cout.flush(); 
	  Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	  ApplyClFactors(tempCl, Nls, lsup, supindex);
	  //cout << "done.\n";
	  // Output:
	  filename=PrintOut(config.reads("REG_CL_PREFIX"), i, j, N1, N2, lls, tempCl, Nls);
	  //cout << ">> Regularized lognormal C(l) written to "+filename<<endl;
	
	  // Temporary memory deallocation:
	  free_vector(tempCl, 0, lastl);
	  free_vector(xi, 0, 2*Nls-1);
	  free_vector(workspace, 0, 2*Nls-1);
	} 
      cout << "done.\n";
      cout << ">> Regularized lognormal C(l) written to prefix "+config.reads("REG_CL_PREFIX")<<endl;
    
    } // End of if(dist==lognormal)
    else warning("ClProcess: regularized C(l)s asked for GAUSSIAN realizations.");
  } // End of computing regularized lognormal Cls.
  

  // Freeing memory: from now on we only need CovByl, means, shifts.
  if (dist==lognormal) {
    cout << "DLT memory deallocation...                                "; cout.flush();
    free_vector(lls, 0, lastl);
    free_vector(LegendreP, 0, 2*Nls*Nls-1);
    free_vector(DLTweights, 0, 4*Nls-1); 
    cout << "done.\n";
  }
  
  // Exit if this is the last output requested:
  if (config.reads("EXIT_AT")=="REG_CL_PREFIX") return 1;

  return 0; // Any return in the middle of this function returns 1.
}
