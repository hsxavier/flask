#include <iostream>
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include <gsl/gsl_linalg.h>     // Cholesky descomposition.
#include <gsl/gsl_randist.h>    // Random numbers.
#include <math.h>               // Log function.
#include <cstdlib>              // For function 'system'.
#include <string>               // For function 'to_string'.

std::ofstream debugfile;        // Send debugging messages to this file.



/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl; using std::string; using std::ofstream; // Basic stuff.
  using namespace ParDef; ParameterList config;                             // Easy configuration file use.
  char message[100];                                                        // Handling warnings and errors.
  std::string filename, temp;
  std::ofstream outfile;                                                    // File for output.
  std::ifstream infile;                                                     // File for input.
  enum simtype {gaussian, lognormal}; simtype dist;                         // For specifying simulation type.
  gsl_matrix *CovMatrix; long CovSize;
  int status, i, j;
  gsl_rng *rnd; double *gaus0, *gaus1;                                      // Random number stuff.
  double *means, *variances, *shifts; long nmeans, nshifts;
  void CorrGauss(double *gaus1, gsl_matrix *L, double *gaus0);
  int GetGaussCov(gsl_matrix *gCovar, gsl_matrix *lnCovar, double *means, double *shifts);
  double Gauss2LNvar(double gvar, double mean, double variance, double shift);
  int getll(const std::string filename);
  std::string getllstr(const std::string filename);
  gsl_set_error_handler_off();                                              // !!! All GSL return messages MUST be checked !!!

  // Opening debug file for dumping information about the program:
  debugfile.open("debug.log");
  if (!debugfile.is_open()) warning("corrlnfields: cannot open debug file.");

  // Loading config file:
  if (argc<=1) { cout << "You must supply a config file." << endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl; 
  
  // Read input data:
  // - Type of simulation
  if (config.reads("DIST")=="LOGNORMAL") dist=lognormal;
  else if (config.reads("DIST")=="GAUSSIAN") dist=gaussian;
  else error("corrlnfields: unknown DIST: "+config.reads("DIST"));
  // - Set random number generator
  rnd = gsl_rng_alloc(gsl_rng_mt19937);
  if (rnd==NULL) error("corrlnfields: gsl_rng_alloc failed!");
  gsl_rng_set(rnd, config.readi("RNDSEED"));    // set random seed
  // - Means and shifts
  means     = LoadList<double>(config.reads("MEANS"), &nmeans);
  if (dist==lognormal)
    shifts  = LoadList<double>(config.reads("SHIFTS"), &nshifts);
  // - Covariance matrices list
  sprintf(message, "ls %s* > corrlnfields.temp", config.reads("COV_PREFIX").c_str());
  system(message);
  infile.open("corrlnfields.temp");
  if (!infile.is_open()) error("corrlnfields: cannot open file corrlnfields.temp");
  
  
  // LOOP over l's:
  while (infile >> filename) {
    cout << "** Working with file '"+filename+"':\n   ";
    CovMatrix = LoadGSLMatrix(filename);
    CovSize   = CovMatrix->size1;
    if (dist==lognormal) {
      variances = vector<double>(0,CovSize-1);
      for (i=0; i<CovSize; i++) variances[i] = CovMatrix->data[i*CovSize+i];  
    }
    
    // Input sanity checks:
    cout << "   Performing sanity checks... "; cout.flush();
    if (CovMatrix->size1 != CovMatrix->size2) error("corrlnfields: Cov. matrix is not square.");
    if (CovSize != nmeans) error("corrlnfields: Cov. matrix and means vector size do not match.");
    if (dist==lognormal) {
    if (CovSize != nshifts) error("corrlnfields: Cov. matrix and shifts vector size do not match.");
    for (i=0; i<CovSize; i++) if(means[i]+shifts[i]<=0) { 
	sprintf(message, "corrlnfields: mean+shift at position %d must be greater than zero.", i); error(message);
      }
    }
    cout << "done.\n";
    
    // If DIST=LOGNORMAL, transform to gaussian variables in place:
    if (dist==lognormal) {
      cout << "   Transforming COV_MATRIX to gaussian... "; cout.flush();
      status=GetGaussCov(CovMatrix, CovMatrix, means, shifts);
      if (status==EDOM) error("corrlnfields: GetGaussCov found bad log arguments.");
      cout << "done.\n";
      // Output gaussian covariance matrix if asked:
      if (config.reads("GCOVOUT_PREFIX")!="0") {
	temp = config.reads("GCOVOUT_PREFIX") + getllstr(filename) + ".dat";
	outfile.open(temp.c_str());
	if (!outfile.is_open()) 
	  warning("corrlnfields: cannot open "+temp+" file.");
	else { 
	  PrintGSLMatrix(CovMatrix, &outfile); outfile.close(); 
	  cout << "   Gaussian covariance matrix written to "+temp << endl;
	}
      }
    }
    
    // Perform a Cholesky decomposition:
    cout << "   Performing a Cholesky decomposition of COV_MATRIX... "; cout.flush();
    status = gsl_linalg_cholesky_decomp(CovMatrix);
    if (status==GSL_EDOM) error("Cholesky decomposition failed: matrix is not positive-definite.");
    cout << "done.\n";
    
    // Generate independent 1sigma random variables:
    cout << "   Generating random variables... "; cout.flush();
    gaus0 = vector<double>(0, CovSize-1);
    gaus1 = vector<double>(0, CovSize-1);
    for (j=0; j<1000; j++) {
      for (i=0; i<CovSize; i++) gaus0[i] = gsl_ran_gaussian(rnd, 1.0); 
      
      // Generate correlated gaussian variables according to CovMatrix:
      CorrGauss(gaus1, CovMatrix, gaus0);
      
      // If DIST=LOGNORMAL, transform variables to lognormal:
      if (dist==lognormal)
	for (i=0; i<CovSize; i++) gaus0[i] = Gauss2LNvar(gaus1[i], means[i], variances[i], shifts[i]);
      else if (dist==gaussian) 
	for (i=0; i<CovSize; i++) gaus0[i] = gaus1[i] + means[i];
      
      sprintf(message,"%d %.8g %.8g %.8g %.8g\n",getll(filename),gaus0[0],gaus0[1],gaus0[2],gaus0[3]);
      debugfile << message;
    }
    cout << "done.\n";
    
    // Free memory allocated inside the loop: avoid memory leakage!
    gsl_matrix_free(CovMatrix);
    free_vector(gaus0,0,CovSize-1);
    free_vector(gaus1,0,CovSize-1);
    if (dist==lognormal) free_vector(variances, 0, CovSize-1);
  } // End of LOOP over l's.

  // End of the program
  infile.close();
  if (dist==lognormal) free_vector(shifts,0,CovSize-1);
  free_vector(means,0,CovSize-1);
  system("rm -f corrlnfields.temp");
  gsl_rng_free(rnd);
  debugfile.close();
  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}



/*** Multiply Lower-triangular matrix L to vector gaus0 and return gaus1 ***/
void CorrGauss(double *gaus1, gsl_matrix *L, double *gaus0) {
  long i, j;

  for(i=0; i<L->size1; i++) {
    gaus1[i]=0;     // L matrix stored as vector in row-major order.
    for(j=0; j<=i; j++) gaus1[i] += L->data[i*L->size1+j] * gaus0[j]; 
  }
}


/*** Transform Cov. matrix of lognormal variables into 
     the Cov. matrix of associated gaussian variables  ***/
int GetGaussCov(gsl_matrix *gCovar, gsl_matrix *lnCovar, double *means, double *shifts) {
  double arg;
  long i, j;
  char message[100]; int status=0; double bad=-666.0;
  
  for(i=0; i<lnCovar->size1; i++)
    for(j=0; j<lnCovar->size2; j++) {
      arg = lnCovar->data[i*(lnCovar->size1)+j]/(means[i]+shifts[i])/(means[j]+shifts[j]) + 1.0;
      // If arg<0, it is an error.
      if (arg<=0) {
	sprintf(message,"GetCovOfGauss: lnCovar(%ld, %ld) leads to bad log argument, gCovar(%ld, %ld) set to %g.",i,j,i,j,bad);
	warning(message);
	status=EDOM;
	gCovar->data[i*(gCovar->size1)+j] = bad;
      }
      // If arg>0, it is valid.
      else gCovar->data[i*(gCovar->size1)+j] = log(arg);  
    }
  return status;
}


/*** Compute a log-normal variable from a gaussian one ***/
double Gauss2LNvar(double gvar, double mean, double variance, double shift) {
  double expsigma2, expmu;
  
  expsigma2 = 1 + variance/pow(mean+shift,2);
  expmu = (mean+shift)/sqrt(expsigma2);

  return expmu*exp(gvar)-shift;
}


/*** Get a number that specify the l of the cov. matrix ***/
int getll(const std::string filename) {
  int i=0, num=0, fileL;
  
  fileL=filename.length();
  // Find a number:
  while (isdigit(filename.c_str()[i])==0) {i++; if(i>=fileL) error("getll: cannot find any number.");}
  // Read the number:
  while (isdigit(filename.c_str()[i])!=0) {num = num*10 + (filename.c_str()[i]-'0'); i++;}
  // Check if there are more numbers in filename:
  while (i<=fileL) {
    if (isdigit(filename.c_str()[i])!=0) error("getll: found more numbers than expected.");
    i++;
  }

  // Return number found:
  return num;
}


/*** Returns getll as string ***/
std::string getllstr(const std::string filename) {
  std::stringstream ss;
  ss << getll(filename);
  return ss.str();
}
