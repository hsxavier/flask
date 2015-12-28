#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling and memory allocation.
#include "s2kit10_naive.hpp"    // For Discrete Legendre Transforms.
#include <gsl/gsl_matrix.h>     // gsl_matrix.
#include <math.h>               // exp, log...
#include "GeneralOutput.hpp"
#include "RegularizeCov.hpp"
#include "flask_aux.hpp"        // For n2fz functions, Maximum, etc..
#include <unistd.h>             // For access function.
#include "FieldsDatabase.hpp"   
#include "fitsfunctions.hpp"    // For ReadHealpixData function used for Healpix pixel window function.
#include "Spline.hpp"           // For applying Healpix window function to arbitrarily spaced C(l)s.


/*** Transform a C(l) to represent the 2D field now smoothed by a Gaussian with variance sigma2 ***/
void ApplyGausWinFunc(double *ClOut, double sigma2, double *l, double *ClIn, int Nls) {
  int i;
  for (i=0; i<Nls; i++)
    ClOut[i] = exp( -l[i]*(l[i]+1)*sigma2 ) * ClIn[i];
} 


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


/*** Export function y(x) for the field combination [i,j] to file ***/
std::string PrintOut(std::string prefix, int i, int j, const FZdatabase & fieldlist, double *x, double *y, int length) {
  int af, az, bf, bz;
  char message[100];
  std::string filename;
  double *wrapper[2];
  std::ofstream outfile;

  wrapper[0] =  x;
  wrapper[1] =  y;

  fieldlist.Index2Name(i, &af, &az); fieldlist.Index2Name(j, &bf, &bz); 
  sprintf(message, "%sf%dz%df%dz%d.dat", prefix.c_str(),af,az,bf,bz);
  filename.assign(message);

  outfile.open(message);
  if (!outfile.is_open()) error("PrintOut: cannot open file "+filename);
  PrintVecs(wrapper, length, 2, &outfile);
  outfile.close();

  return filename;
}


/*** Wraps all processing of the input Cls up to the gaussian covariance matrices for each l ***/
int ClProcess(gsl_matrix ***CovBylAddr, int *NlsOut, const FZdatabase & fieldlist, const ParameterList & config) {
  using namespace definitions;                          // Global definitions.
  using std::cout; using std::endl;                     // Basic stuff.
  simtype dist;                                         // For specifying simulation type.
  char message[200];                                    // Writing outputs with sprintf.
  std::ofstream outfile;                                // File for output.
  FILE* stream; int NinputCls; std::string *filelist;   // To list input Cls.
  int i, j, k, l, m, status, Nfields, Nf, Nz, Nls, lmin, lmax;
  std::string filename, ExitAt, prefix;
  bool *fnzSet;
  gsl_matrix **CovByl;
  double temp, badcorrfrac;

  // Getting general information:
  Nfields = fieldlist.Nfields();
  if (config.reads("DIST")=="LOGNORMAL") dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN") dist=gaussian;
  ExitAt  = config.reads("EXIT_AT");
  lmax    = config.readi("LRANGE", 1);
  lmin    = config.readi("LRANGE", 0);
  if (lmax<lmin) error("ClProcess: LRANGE set in the wrong order.");

  /********************************************/
  /*** PART 1: Load C(l)s and organize them ***/
  /********************************************/
  const int HWMAXL = 10000000; int lastl = HWMAXL;
  int af, az, bf, bz, Nlinput, **fnz, **NentMat;
  long *Nentries, ncols;
  double ***ll, ***Cov, *wrapper[2];
  bool **IsSet;

  // Get list of the necessary C(l) files:
  prefix    = config.reads("CL_PREFIX"); 
  NinputCls = Nfields*Nfields;
  filelist  = vector<std::string>(0,NinputCls-1);
  // LOOP over all C(l)s:
  for (k=0; k<NinputCls; k++) {
    i = k/Nfields;      j = k%Nfields;
    fieldlist.Index2Name(i, &af, &az);
    fieldlist.Index2Name(j, &bf, &bz);
    sprintf(message, "%sf%dz%df%dz%d.dat", prefix.c_str(), af, az, bf, bz);
    if(access(message, R_OK) == 0) filelist[k].assign(message);
  }
  
  // Find out the number of lines and columns in C(l) files:  
  Nlinput  = 0;
  Nentries = vector<long>(0,NinputCls-1);
  for (k=0; k<NinputCls; k++) {
    if (filelist[k].size()>0) {
      CountEntries(filelist[k], &(Nentries[k]), &ncols); // Get number of Nls.
      if (ncols!=2) error("ClProcess: wrong number of columns in file "+filelist[k]);
      if (Nentries[k]>Nlinput) Nlinput=Nentries[k];          // Record maximum number of ls.
    }
    else Nentries[k]=0;
  }
  
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
  for (k=0; k<NinputCls; k++) if (filelist[k].size()>0) {
      i = k/Nfields;      j = k%Nfields;
      fieldlist.Index2Name(i, &af, &az);
      fieldlist.Index2Name(j, &bf, &bz);      
      cout << filelist[k] << " goes to ["<<i<<", "<<j<<"]" << endl;
      // Record the order of the fields in CovMatrix:  
      if (fnzSet[i]==0) { fnz[i][0] = af; fnz[i][1] = az; fnzSet[i] = 1; }
      else if (fnz[i][0] != af || fnz[i][1] != az) error("ClProcess: field order in CovMatrix is messed up!"); 
      if (fnzSet[j]==0) { fnz[j][0] = bf; fnz[j][1] = bz; fnzSet[j] = 1; }
      else if (fnz[j][0] != bf || fnz[j][1] != bz) error("ClProcess: field order in CovMatrix is messed up!");
      // Import data:
      wrapper[0] = &(ll[i][j][0]);
      wrapper[1] = &(Cov[i][j][0]);
      ImportVecs(wrapper, Nentries[k], 2, filelist[k].c_str());
      NentMat[i][j] = Nentries[k];
      IsSet[i][j]=1; 
    };
  if (config.readi("ALLOW_MISS_CL")==1) cout << "ALLOW_MISS_CL=1: will set totally missing Cl's to zero.\n";
  free_vector(Nentries, 0, NinputCls-1);
  free_vector(filelist, 0, NinputCls-1);

  // Check if every field was assigned a position in the CovMatrix:
  for (i=0; i<Nfields; i++) if (fnzSet[i]==0) error("ClProcess: some position in CovMatrix is unclaimed.");
  free_vector(fnzSet, 0, Nfields-1);
  // If positions are OK and output required, print them out:
  if (config.reads("FLIST_OUT")!="0") {
    outfile.open(config.reads("FLIST_OUT").c_str());
    if (!outfile.is_open()) error("ClProcess: cannot open FLIST_OUT file.");
    PrintTable(fnz, Nfields, 2, &outfile);
    outfile.close();
    cout << ">> Written field list to "+config.reads("FLIST_OUT")<<endl;
  }
  free_matrix(fnz, 0, Nfields-1, 0,1);
  // Exit if this is the last output requested:
  if (ExitAt=="FLIST_OUT") return 1;
 
  
  /*************************************************************/
  /*** PART 1.5: Apply various window functions if requested ***/
  /*************************************************************/
  double WinFuncVar, *pixwin, *pixell;
  Spline pixSpline;
  

  // Gaussian beam:
  WinFuncVar = config.readd("WINFUNC_SIGMA");            // WINFUNC_SIGMA will be transformed to radians and squared below. 
  if (WinFuncVar > 0.0) {
    Announce("Applying Gaussian window function to C(l)s... ");
    WinFuncVar = WinFuncVar/60.0*M_PI/180.0;             // From arcmin to radians.
    WinFuncVar = WinFuncVar*WinFuncVar;                  // From std. dev. to variance.
    // LOOP over existing C(l)s:
#pragma omp parallel for schedule(dynamic) private(i, j)
    for (k=0; k<Nfields*Nfields; k++) {
      i=k/Nfields;  j=k%Nfields;
      if (IsSet[i][j]==1) {
	// In-place C(l) change due to Gaussian window function:
	ApplyGausWinFunc(Cov[i][j], WinFuncVar, ll[i][j], Cov[i][j], NentMat[i][j]);
      }
    } // End over LOOP over existing C(l)s.
    Announce();
  } // End of IF Smoothing requested.


  // Healpix pixel window function:      
  if (config.readi("APPLY_PIXWIN")==1) {
    Announce("Applying Healpix pixel window function to C(l)s... ");
    // Prepare spline of the window function [input C(l)s might be at random ell]:
    m      = config.readi("NSIDE");
    pixell = vector<double>(0, 4*m);
    pixwin = vector<double>(0, 4*m);
    status = ReadHealpixData(1, config, pixwin, 2);
    if (status!=0) error("ClProcess: cannot read Healpix pixel window FITS.");
    for (i=0; i<=4*m; i++) {
      pixell[i] = (double)i;
      pixwin[i] = pixwin[i]*pixwin[i];
    }
    pixSpline.init(pixell, pixwin, 4*m+1);
    free_vector(pixell, 0, 4*m);
    free_vector(pixwin, 0, 4*m);
    // LOOP over existing C(l)s:
#pragma omp parallel for schedule(dynamic) private(i, j, l)
    for (k=0; k<Nfields*Nfields; k++) {
      i=k/Nfields;  j=k%Nfields;
      if (IsSet[i][j]==1) {
	// In-place C(l) change due to pixel window function:
	if(ll[i][j][NentMat[i][j]-1] > 4*m) warning("ClProcess: input C(l) overshoot Healpix ppixel window function.");
	for(l=0; l<NentMat[i][j]; l++) Cov[i][j][l] = pixSpline(ll[i][j][l])*Cov[i][j][l];
      }
    } // End over LOOP over existing C(l)s.
    Announce();
  } // End of IF Smoothing requested.


  // Print C(l)s to files if requested:
  filename = config.reads("SMOOTH_CL_PREFIX");
  if (filename!="0") {
    for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]==1) {
	  PrintOut(filename, i, j, fieldlist, ll[i][j], Cov[i][j], NentMat[i][j]); 	  
	}
    cout << ">> Smoothed C(l)s written to prefix "+config.reads("SMOOTH_CL_PREFIX")<<endl;
  }
  if (ExitAt=="SMOOTH_CL_PREFIX") return 1;


  /*** Continue organizing C(l)s ***/

  // Look for the maximum l value described by all C(l)s:
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]==1) {
	if (ll[i][j][NentMat[i][j]-1]>HWMAXL) error ("ClProcess: too high l in C(l)s: increase HWMAXL.");
	if (ll[i][j][NentMat[i][j]-1]<lastl) lastl = (int)ll[i][j][NentMat[i][j]-1];
      }
  Nls=lastl+1; // l=0 is needed for DLT. Nls is known as 'bandwidth' (bw) in s2kit 1.0 code.
  (*NlsOut)=Nls;

  
  // Allocate gsl_matrices that will receive covariance matrices for each l.
  Announce("Allocating data-cube needed for Cholesky decomposition... ");
  CovByl = GSLMatrixArray(Nls, Nfields, Nfields);
  (*CovBylAddr) = CovByl;
  Announce();
 

  /*****************************************************************/
  /*** PART 2: Compute auxiliary gaussian C(l)s if LOGNORMAL     ***/
  /*****************************************************************/
  double *tempCl, *LegendreP, *workspace, *xi, lsup, supindex, *theta, *DLTweights, *lls;

  if (dist==lognormal) {
    cout << "LOGNORMAL realizations: will compute auxiliary gaussian C(l)s:\n";
    // Loads necessary memory:
    Announce("Allocating memory for DLT... ");
    lls        = vector<double>(0, lastl);
    workspace  = vector<double>(0, 16*Nls-1);
    LegendreP  = vector<double>(0, 2*Nls*Nls-1);
    DLTweights = vector<double>(0, 4*Nls-1);
    // Initialize vectors:
    for (i=0; i<=lastl; i++) lls[i]=(double)i;
    Announce();// angle theta is only necessary for output:
    if (config.reads("XIOUT_PREFIX")!="0" || config.reads("GXIOUT_PREFIX")!="0") {
      Announce("Generating table of sampling angles... ");
      theta    = vector<double>(0, 2*Nls-1);
      ArcCosEvalPts(2*Nls, theta);
      for (i=0; i<2*Nls; i++) theta[i] = theta[i]*180.0/M_PI;
      Announce();
    } 
    
    // Loads C(l) exponential suppression:
    lsup     = config.readd("SUPPRESS_L");
    supindex = config.readd("SUP_INDEX"); 
    // Load s2kit 1.0 Legendre Polynomials:
    Announce("Generating table of Legendre polynomials... ");
    PmlTableGen(Nls, 0, LegendreP, workspace);
    free_vector(workspace, 0, 16*Nls-1);
    Announce();
    // Compute s2kit 1.0 Discrete Legendre Transform weights:
    Announce("Calculating forward DLT weights... ");
    makeweights(Nls, DLTweights);
    Announce();
  }

  // LOOP over all C(l)s already set.
  if (dist==lognormal) Announce("Transforming C(l)s for the auxiliary Gaussian ones... ");
  else Announce("Interpolating C(l)s for all l's... ");
#pragma omp parallel for schedule(dynamic) private(tempCl, xi, workspace, filename, l, i, j)
  for (k=0; k<Nfields*Nfields; k++) {
    i=k/Nfields;  j=k%Nfields;
    if (IsSet[i][j]==1) {
      
      // Temporary memory allocation:
      tempCl    = vector<double>(0, lastl);
      xi        = vector<double>(0, 2*Nls-1);
      workspace = vector<double>(0, 2*Nls-1);

      // Interpolate C(l) for every l; input C(l) might not be like that:
      GetAllLs(ll[i][j], Cov[i][j], NentMat[i][j], tempCl, lastl, config.readi("EXTRAP_DIPOLE"));
	
      if (dist==lognormal) {              /** LOGNORMAL ONLY **/
	
	// Compute correlation function Xi(theta):
	ModCl4DLT(tempCl, lastl, lsup, supindex);
	Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	if (config.reads("XIOUT_PREFIX")!="0") { // Write it out if requested:
	  filename=PrintOut(config.reads("XIOUT_PREFIX"), i, j, fieldlist, theta, xi, 2*Nls);
	}

	// Transform Xi(theta) to auxiliary gaussian Xi(theta):
	status=GetGaussCorr(xi, xi, 2*Nls, fieldlist.mean(i), fieldlist.shift(i), fieldlist.mean(j), fieldlist.shift(j));
	if (status==EDOM) error("ClProcess: GetGaussCorr found bad log arguments.");
	if (i==j && xi[0]<0) warning("ClProcess: auxiliary field variance is negative.");
	if (config.reads("GXIOUT_PREFIX")!="0") { // Write it out if requested:
	  filename=PrintOut(config.reads("GXIOUT_PREFIX"), i, j, fieldlist, theta, xi, 2*Nls);
	}

	// Transform Xi(theta) back to C(l):
	Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	ApplyClFactors(tempCl, Nls);
	if (config.reads("GCLOUT_PREFIX")!="0") { // Write it out if requested:
	  filename=PrintOut(config.reads("GCLOUT_PREFIX"), i, j, fieldlist, lls, tempCl, Nls);
	}	  
      }                                 /** END OF LOGNORMAL ONLY **/ 
	
      // Save auxiliary C(l):
      for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j]=tempCl[l];
	
      // Temporary memory deallocation:
      free_vector(tempCl, 0, lastl);
      free_vector(xi, 0, 2*Nls-1);
      free_vector(workspace, 0, 2*Nls-1);
    } // End of IF C(l)[i,j] is set.
  } // End of LOOP over C(l)[i,j] that were set.
  Announce();
  
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
  if (ExitAt=="XIOUT_PREFIX"  || 
      ExitAt=="GXIOUT_PREFIX" || 
      ExitAt=="GCLOUT_PREFIX") return 1;
    
  // Set Cov(l)[i,j] = Cov(l)[j,i]
  Announce("Set remaining cov. matrices elements based on symmetry... ");
  k = config.readi("ALLOW_MISS_CL");
  if (k!=1 && k!=0) error("ClProcess: unknown option for ALLOW_MISS_CL.");
  for(i=0; i<Nfields; i++)
    for(j=0; j<Nfields; j++)
      // Look for empty entries in Cov:
      if (IsSet[i][j]==0) {
	// If transpose is empty too:
	if (IsSet[j][i]==0) {
	  // Set transpose to zero if this is allowed:
	  if (k==1) { for (l=0; l<Nls; l++) CovByl[l]->data[j*Nfields+i] = 0.0; IsSet[j][i]=1; }
	  // If not allowed, return error:
	  else { sprintf(message,"ClProcess: [%d, %d] could not be set because [%d, %d] was not set.",i,j,j,i); error(message); }
	}
	for (l=0; l<Nls; l++) CovByl[l]->data[i*Nfields+j] = CovByl[l]->data[j*Nfields+i];
	IsSet[i][j] = 1;
      }
  for(i=0; i<Nfields; i++) for(j=0; j<Nfields; j++) if (IsSet[i][j]!=1) {
	sprintf(message,"ClProcess: [%d, %d] was not set.",i,j);
	error(message);
      }
  Announce();
  free_matrix(IsSet, 0, Nfields-1, 0, Nfields-1);

  // Output covariance matrices for each l if requested:
  GeneralOutput(CovByl, config, "COVL_PREFIX", 0);
  if (config.reads("COVL_PREFIX")!="0") 
    cout << ">> Cov. matrices written to prefix "+config.reads("COVL_PREFIX")<<endl;
  // Exit if this is the last output requested:
  if (ExitAt=="COVL_PREFIX") return 1;
  
  // Verify basic properties of auxiliary cov. matrices:
  Announce("Verifying aux. Cov. matrices properties... ");
  badcorrfrac = config.readd("BADCORR_FRAC");
  for (l=lmin; l<=lmax; l++) // Skipping l=0 since matrix should be zero.
    for (i=0; i<Nfields; i++) {
      // Verify that diagonal elements are positive:
      if (CovByl[l]->data[i*Nfields+i]<=0.0) {
	sprintf(message, "ClProcess: Cov. matrix (l=%d) element [%d, %d] is negative or zero", l, i, i);
	warning(message);
      }
      for (j=i+1; j<Nfields; j++) {
	// Correlations c should be limited to -1<c<1.
	temp = CovByl[l]->data[i*Nfields+j]/sqrt(CovByl[l]->data[i*Nfields+i]*CovByl[l]->data[j*Nfields+j]);
	if (temp>1.0 || temp<-1.0) {
	  // Try increasing variances if correlation is absurd:
	  cout << "  Aux. Cov. matrix (l="<<l<<") element ["<<i<<", "<<j<<"] results in correlation "<<temp
	       <<". Fudging variances with BADCORR_FRAC...\n";
	  CovByl[l]->data[i*Nfields+i] *= (1.0+badcorrfrac);
	  CovByl[l]->data[j*Nfields+j] *= (1.0+badcorrfrac);
	  temp = CovByl[l]->data[i*Nfields+j]/sqrt(CovByl[l]->data[i*Nfields+i]*CovByl[l]->data[j*Nfields+j]);
	  if (temp>1.0 || temp<-1.0) warning("ClProcess: BADCORR_FRAC could not solve the issue.");
	}
      }
    }      
  Announce();


  /****************************************************************************/
  /*** PART 3: Regularize (make them positive definite) covariance matrices ***/
  /****************************************************************************/
  gsl_matrix *gslm;
  double *MaxChange, MMax;
  int lMMax, lstart, lend, FailReg=0, NCls;

  // If producing the regularized lognormal C(l)s, all ls must be regularized (skip l=0 because is set to zero).
  // Else, only the cov. matrices for requested ls are regularized.
  // Note that a very high exponential suppression of C(l)s make difficult to regularize matrices.
  if (dist==lognormal && config.reads("REG_CL_PREFIX")!="0") { lstart = 1;  lend = Nls-1; }
  else { lstart = lmin;  lend = lmax; }
  MaxChange = vector<double>(lstart, lend);
  
  Announce("Regularizing cov. matrices... ");
  
  #pragma omp parallel for schedule(dynamic) private(gslm, filename)  
  for (l=lstart; l<=lend; l++) {

    // Check pos. defness, regularize if necessary, keep track of changes:
    gslm = gsl_matrix_alloc(Nfields, Nfields);
    gsl_matrix_memcpy(gslm, CovByl[l]);
    status       = RegularizeCov(CovByl[l], config);
    MaxChange[l] =   MaxFracDiff(CovByl[l], gslm);
    if (status==9) { 
      sprintf(message, "ClProcess: RegularizeCov for l=%d reached REG_MAXSTEPS with Max. change of %g.",l,MaxChange[l]); 
      warning(message);
      FailReg=1;
    }
    gsl_matrix_free(gslm);
    // Output regularized matrix if requested:
    if (config.reads("REG_COVL_PREFIX")!="0") {
      filename=config.reads("REG_COVL_PREFIX")+"l"+ZeroPad(l,lend)+".dat";
      GeneralOutput(CovByl[l], filename, 0);
    }
  }
  Announce();
  if (FailReg==1) error("ClProcess: failed to regularize covariance matrices.");
  
  // Dump changes in cov. matrices to the screen:
  MMax  = 0.0; lMMax = 0;
  for (l=lmin; l<=lmax; l++) if (MaxChange[l]>MMax) {MMax = MaxChange[l]; lMMax = l;}
  cout << "Max. frac. change for "<<lmin<<"<=l<="<<lmax<<" at l="<<lMMax<<": "<<MMax<<endl;  
  free_vector(MaxChange, lstart, lend);
  // Output regularized matrices if requested:
  if (config.reads("REG_COVL_PREFIX")!="0") 
    cout << ">> Regularized cov. matrices written to prefix "+config.reads("REG_COVL_PREFIX")<<endl;
  // Exit if this is the last output requested:
  if (ExitAt=="REG_COVL_PREFIX") return 1;

  /***********************************************************/
  /*** PART 4: Obtain regularized input Cls if requested   ***/
  /***********************************************************/

  if (config.reads("REG_CL_PREFIX")!="0") {
    if(dist==lognormal) {
      Announce("Computing regularized lognormal Cls... ");
            
      // LOOP over fields:
      NCls = (Nfields*(Nfields+1))/2;
#pragma omp parallel for schedule(dynamic) private(tempCl, xi, workspace, filename, l, m, i, j)
      for (k=0; k<NCls; k++) {
	l = (int)((sqrt(8.0*(NCls-1-k)+1.0)-1.0)/2.0);
	m = NCls-1-k-(l*(l+1))/2;
	i = Nfields-1-l;
	j = Nfields-1-m;

	// Temporary memory allocation:
	tempCl    = vector<double>(0, lastl);
	xi        = vector<double>(0, 2*Nls-1);
	workspace = vector<double>(0, 2*Nls-1);

	// Copy the Cl to a vector:
	for (l=0; l<Nls; l++) tempCl[l] = CovByl[l]->data[i*Nfields+j]; // tudo certo.
	// Compute correlation function Xi(theta):
	ModCl4DLT(tempCl, lastl, -1, -1); // Suppression not needed (it was already suppressed).
	Naive_SynthesizeX(tempCl, Nls, 0, xi, LegendreP);
	// Get Xi(theta) for lognormal variables:
	GetLNCorr(xi, xi, 2*Nls, fieldlist.mean(i), fieldlist.shift(i), fieldlist.mean(j), fieldlist.shift(j));
	// Compute the Cls:
	Naive_AnalysisX(xi, Nls, 0, DLTweights, tempCl, LegendreP, workspace);
	ApplyClFactors(tempCl, Nls, lsup, supindex);
	// Output:
	filename=PrintOut(config.reads("REG_CL_PREFIX"), i, j, fieldlist, lls, tempCl, Nls);
	
	// Temporary memory deallocation:
	free_vector(tempCl, 0, lastl);
	free_vector(xi, 0, 2*Nls-1);
	free_vector(workspace, 0, 2*Nls-1);
      } 
      Announce();
      cout << ">> Regularized lognormal C(l) written to prefix "+config.reads("REG_CL_PREFIX")<<endl;
    
    } // End of if(dist==lognormal)
    else warning("ClProcess: regularized C(l)s asked for GAUSSIAN realizations.");
  } // End of computing regularized lognormal Cls.
  

  // Freeing memory: from now on we only need CovByl, means, shifts.
  if (dist==lognormal) {
    Announce("DLT memory deallocation... ");
    free_vector(lls, 0, lastl);
    free_vector(LegendreP, 0, 2*Nls*Nls-1);
    free_vector(DLTweights, 0, 4*Nls-1); 
    Announce();
  }
  
  // Exit if this is the last output requested:
  if (ExitAt=="REG_CL_PREFIX") return 1;

  return 0; // Any return in the middle of this function returns 1.
}
