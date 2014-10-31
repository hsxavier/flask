#include <iostream>
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include <gsl/gsl_linalg.h>     // Cholesky descomposition.
#include <gsl/gsl_randist.h>    // Random numbers.

std::ofstream debugfile;        // Send debugging messages to this file.



/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl; using std::string; using std::ofstream; // Basic stuff.
  using namespace ParDef; ParameterList config;                             // Easy configuration file use.  
  gsl_matrix *CovMatrix; long CovSize;
  int status, i, j;
  gsl_rng *rnd; double *gaus0, *gaus1;                                      // Random number stuff.
  void CorrGauss(double *gaus1, gsl_matrix *L, double *gaus0);
  gsl_set_error_handler_off();                                              // !!! All GSL return messages MUST be checked !!!

  /*** Opening debug file for dumping information about the program ***/
  debugfile.open("debug.log");
  if (!debugfile.is_open()) error("main: cannot open debug file.");

  /*** Loading config file ***/
  if (argc<=1) { cout << "You must supply a config file." << endl; return 0;}
  config.load(argv[1]);
  cout << endl;
  cout << "-- Configuration setup:\n";
  cout << "   File: "<<argv[1]<<endl;
  config.lineload(argc, argv);
  config.show();
  cout << endl; 

  // Read Covariance matrix:
  CovMatrix = LoadGSLMatrix(config.reads("COV_MATRIX"));
  CovSize   = CovMatrix->size1;
  // Perform a Cholesky decomposition:
  status = gsl_linalg_cholesky_decomp(CovMatrix);
  if (status==GSL_EDOM) error("Cholesky decomposition failed: matrix is not positive-definite.");
  
  //PrintGSLMatrix(CovMatrix);
  
  // Generate independent 1sigma random variables:
  gaus0 = vector<double>(0, CovSize-1);
  gaus1 = vector<double>(0, CovSize-1);
  rnd = gsl_rng_alloc(gsl_rng_mt19937);
  if (rnd==NULL) error("main: gsl_rng_alloc failed!");
  gsl_rng_set(rnd, config.readi("RNDSEED"));    // set random seed
  for (j=0; j<1000000; j++) {
    for (i=0; i<CovSize; i++) gaus0[i] = gsl_ran_gaussian(rnd, 1.0); 
    //Generate correlated gaussian variables according to CovMatrix:
    CorrGauss(gaus1, CovMatrix, gaus0);
    debugfile << gaus1[0] << " " << gaus1[1] << " " << gaus1[2] << endl;
  }
  
  

  /*** End of the program ***/
  free_vector(gaus0,0,CovSize-1);
  free_vector(gaus1,0,CovSize-1);
  gsl_rng_free(rnd);
  gsl_matrix_free(CovMatrix);
  debugfile.close(); // Close debug file.
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
