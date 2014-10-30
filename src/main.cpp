#include <iostream>
#include "ParameterList.hpp"    // Configuration and input system.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "gsl_aux.hpp"          // Using and reading GSL matrices.
#include <gsl/gsl_linalg.h>     // Cholesky descomposition.

std::ofstream debugfile;



/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl; using std::string; using std::ofstream; // Basic stuff.
  using namespace ParDef; ParameterList config;                             // Easy configuration file use.  
  gsl_matrix *CovMatrix; long CovSize;
  int status;
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
  // Perform a Cholesky decomposition:
  status = gsl_linalg_cholesky_decomp(CovMatrix);
  if (status==GSL_EDOM) error("Cholesky decomposition failed: matrix is not positive-definite.");
  
  PrintGSLMatrix(CovMatrix, &debugfile);
  

  /*** End of the program ***/
  gsl_matrix_free(CovMatrix);
  debugfile.close(); // Close debug file.
  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}

