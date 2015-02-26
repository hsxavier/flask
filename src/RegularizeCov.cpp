#include "RegularizeCov.hpp"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include "Utilities.hpp"
#include "gsl_aux.hpp" // For PrintGSLMatrix.
#include <cmath>       // For sqrt.

// Normalize elements of a SYMMETRIC matrix by the quadratic sum of its independent elements:
void NormMatrix(gsl_matrix *A) {
  int i, j;
  double norm=0;
  
  // Compute the norm through a quadratic sum of the upper triangular part of matrix:
  for (i=0; i<A->size1; i++)
    for (j=i; j<A->size2; j++) norm += A->data[i*A->size1+j]*A->data[i*A->size1+j];
  norm = sqrt(norm);
  
  // Normalize the whole matrix:
  for (i=0; i<A->size1; i++)
    for (j=0; j<A->size2; j++)
      A->data[i*A->size1+j] = A->data[i*A->size1+j] / norm;
  
}


// Returns an evaluation of the change from vec0 to vec1.
// Current implementation returns the total difference between components of vec1 and vec0 that 
// did not stayed positive.
double GradeVecChange(gsl_vector *vec1, gsl_vector *vec0) {
  int i;
  double sum=0;
  if (vec1->size != vec0->size) error("GradeVecChange: encountered vectors with different sizes.");
  for (i=0; i<vec1->size; i++) 
    if(vec1->data[i]<0 || vec0->data[i]<0) 
      sum += vec1->data[i] - vec0->data[i];
  return sum;
}


// Computes the combination of changes in 'input' matrix elements that maximizes the change in eigenvalues: 
void GetDirection(gsl_matrix *direction, double step, gsl_matrix *input) {
  int status, i, j;
  gsl_matrix *onedir;
  gsl_vector *eval0, *eval1;
  gsl_eigen_symm_workspace *workspace;
  
  // Allocate memory:
  workspace = gsl_eigen_symm_alloc(input->size1);
  eval0     = gsl_vector_alloc(input->size1);
  eval1     = gsl_vector_alloc(input->size1);
  onedir    = gsl_matrix_alloc(input->size1, input->size2);
  
  // Compute the original eigenvalues of input:
  status    = gsl_matrix_memcpy(onedir, input);
  status    = gsl_eigen_symm(onedir, eval0, workspace); // This destroys onedir!
  if (status != GSL_SUCCESS) error("GetDirection: cannot compute eigenvalues for input matrix.");
  
  // Select one independent element of matrix, change it and record the effect over negative eigenvalues:
  for (i=0; i<input->size1; i++)
    for (j=i; j<input->size2; j++) {
      // Get copy of original matrix:
      status = gsl_matrix_memcpy(onedir, input);
      // Shift the element:
      onedir->data[i*onedir->size1+j] = onedir->data[i*onedir->size1+j] * (1.0 + step);
      onedir->data[j*onedir->size1+i] = onedir->data[i*onedir->size1+j];      
      // Compute eigenvalue:
      status = gsl_eigen_symm(onedir, eval1, workspace); // This destroys onedir!
      if (status != GSL_SUCCESS) error("GetDirection: cannot compute eigenvalues for onedir matrix.");
      // Set direction element according to the change in negative eigenvalues:
      direction->data[i*direction->size1+j] = GradeVecChange(eval1, eval0);
      direction->data[j*direction->size1+i] = direction->data[i*direction->size1+j];
    }
  // Normalize direction:  
  NormMatrix(direction);

  // Free memory:
  gsl_eigen_symm_free(workspace);
  gsl_vector_free(eval0);
  gsl_vector_free(eval1);
  gsl_matrix_free(onedir);
}


// Turns a symmetric matrix A into a positive definite matrix with the least possible change: 
void RegularizeCov(gsl_matrix * A, const ParameterList & config) {
  int method, status, i, j, k;
  double NewEval;
  const int none=0, frobenius=1, stepper=2;
  gsl_vector *eval;
  gsl_eigen_symm_workspace *workspace;
  gsl_matrix *testmatrix;

  method  = config.readi("REGULARIZE_METHOD");
  NewEval = config.readd("NEW_EVAL"); 
  if (method>2 || method<0) error("RegularizeCov: unkown method.");

  // Run regularization if asked for:
  if (method != none) {
    testmatrix = gsl_matrix_alloc (A->size1, A->size2);
    status     = gsl_matrix_memcpy(testmatrix, A);
    
    // Test if matrix is positive definite and only proceed with regularization if otherwise:
    workspace = gsl_eigen_symm_alloc(testmatrix->size1);
    eval      = gsl_vector_alloc(testmatrix->size1);
    status    = gsl_eigen_symm(testmatrix, eval, workspace); // This DESTROYS input matrix !!!
    gsl_eigen_symm_free(workspace);
    gsl_matrix_free(testmatrix);
    if (status != GSL_SUCCESS) error("RegularizeCov: cannot get matrix eigenvalues.");
    // Look for negative eigenvalue:
    status=0;
    for (i=0; i<eval->size && status!=1; i++) if (gsl_vector_get(eval,i)<0) status=1;
    
    // Regularization is necessary and will be performed:
    if (status==1) {
      
      // Verify if input matrix has the necessary properties:
      if (A->size1 != A->size2) error("RegularizeCov: matrix is not square.");
      for(i=0; i<A->size1; i++) 
	for (j=i+1; j<A->size2; j++) 
	  if (A->data[i*A->size1+j] != A->data[j*A->size2+i]) error("RegularizeCov: matrix is not symmetric.");


      // METHOD: Use positive-semidefinite approximant matrix that minimizes the elements 
      //         differences quadratic sum (Frobenius norm):
      //         This corresponds to setting the negative eigenvalues to zero.
      if (method==frobenius) {
	gsl_eigen_symmv_workspace *workspaceV;
	gsl_matrix *evec;

	std::cout << "   Regularizing matrix by minimizing the elements quadratic sum... "; std::cout.flush();
	// Allocate memory (vector for eigenvalues was already allocated):
	workspaceV = gsl_eigen_symmv_alloc(A->size1);
	evec       = gsl_matrix_alloc (A->size1, A->size1);
	// Get eigenvectors and eigenvalues:
	status     = gsl_eigen_symmv(A, eval, evec, workspaceV);
	gsl_eigen_symmv_free(workspaceV);
	if (status != GSL_SUCCESS) error("RegularizeCov: cannot get matrix eigenvectors and eigenvalues.");
	// Clear negative eigenvalues:
	for (i=0; i<eval->size; i++) if (gsl_vector_get(eval,i)<0) gsl_vector_set(eval, i, NewEval);
	// Compute regularized matrix:
	for(i=0; i<A->size1; i++) 
	  for (j=0; j<A->size2; j++) {
	    gsl_matrix_set(A, i,j, 0.0);
	    for (k=0; k<evec->size2; k++)
	      A->data[i*A->size1+j] += gsl_matrix_get(evec,i,k) * gsl_vector_get(eval,k) * gsl_matrix_get(evec,j,k);
	  }
	gsl_matrix_free(evec);
	std::cout << "done.\n";
      } // End of minimization of frobenius norm method.
      
      
      // METHOD: Make successful approximations trying to get a positive-definite matrix 
      //         with the minimum fractional change in the original matrix:
      if (method==stepper) {
	double step;
	bool posdef;
	gsl_matrix *direction;

	std::cout << "   Regularizing matrix by successive approximations... "; std::cout.flush();
	step = config.readd("REGULARIZE_STEP");
	
	// Allocate memory:
	workspace  = gsl_eigen_symm_alloc(A->size1);
	direction  = gsl_matrix_alloc(A->size1, A->size2);
	testmatrix = gsl_matrix_alloc (A->size1, A->size2);

	// While distorted matrix is not positive-definite, keep on distorting:
	posdef = 0;
	while (posdef==0) {
	  // Find direction with most change in eigenvalues:
	  GetDirection(direction, step, A);
	  // Distorts matrix 'A' in that direction and magnitude 'step':
	  for (i=0; i<A->size1; i++)
	    for (j=i; j<A->size2; j++) {
	      A->data[i*A->size1+j] = A->data[i*A->size1+j] * (1.0 + step*direction->data[i*A->size1+j]);
	      A->data[j*A->size1+i] = A->data[i*A->size1+j]; // A is symmetric.
	    }
	  // Verify if all eigenvalues have become positive:
	  status = gsl_matrix_memcpy(testmatrix, A);
	  status = gsl_eigen_symm(testmatrix, eval, workspace);
	  if (status != GSL_SUCCESS) error("RegularizeCov: cannot get eigenvalues for distorted matrix.");
	  for (i=0; i<eval->size && eval->data[i]>0; i++);
	  // Exit if that is true.
	  if(i >= eval->size) posdef = 1;

	} // End of LOOP while eigenvalues are negative.
	gsl_matrix_free(testmatrix);
	gsl_eigen_symm_free(workspace);
	gsl_matrix_free(direction);
	std::cout << "done.\n";
      } // End of minimum fractional difference method.


    } // End of "regularization is necessary" block.  
    gsl_vector_free(eval);
  } // End of "regularization YES" block.
}


// Compute the matrix eigenvalues without destroying the matrix (by making a copy first):
void GetEigenValues(gsl_vector *eval, gsl_matrix *A, gsl_eigen_symm_workspace *w) {
  int status;
  gsl_matrix *M;
  
  M      = gsl_matrix_alloc(A->size1, A->size2);
  status = gsl_matrix_memcpy(M, A);
  status = gsl_eigen_symm(M, eval, w);
  if (status != GSL_SUCCESS) error("GetEigenValues: cannot get matrix eigenvalues.");
  gsl_matrix_free(M);
}
