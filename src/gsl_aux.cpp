#include <iostream>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include "Utilities.hpp"
#include "gsl_aux.hpp"



/*** Allocate an array of gsl_matrices ***/
gsl_matrix **GSLMatrixArray(int Nmatrices, int Nrows, int Ncols) {
  gsl_matrix **array;
  int i;
  char message[100];
  
  array = (gsl_matrix**) malloc(sizeof(gsl_matrix*)*Nmatrices);
  if (array==NULL) error("GSLMatrixArray: failed to allocate array of gsl_matrix pointers.");
  for(i=0; i<Nmatrices; i++) {
    array[i] = gsl_matrix_alloc(Nrows, Ncols);
    if (array[i]==NULL) {
      sprintf(message,"GSLMatrixArray: failed to allocate i=%d gsl_matrix.", i);
      error(message);
    }
  }
  return array;
}


/*** Free memory allocated for an array of gsl_matrices ***/
void free_GSLMatrixArray(gsl_matrix **array, int Nmatrices) {
  int i;
  
  for(i=0; i<Nmatrices; i++) gsl_matrix_free(array[i]);
  free(array);
}


/*** Import a matrix from a file to a GSLMatrix already allocated ***/
void LoadGSLMatrix(std::string filename, gsl_matrix *matrix) {
  FILE *file;
  int status;

  // Open file:
  if ((file=fopen(filename.c_str(), "r"))==NULL) error("LoadGSLMatrix: cannot open file "+filename);
  
  // Loading values to table:
  status=gsl_matrix_fscanf(file, matrix);
  if(status==GSL_EFAILED) error("LoadGSLMatrix: gsl_matrix_fscanf failed!");
  
  fclose(file);  
}


/*** Import a table from file to gsl_matrix format ***/
gsl_matrix *LoadGSLMatrix(std::string filename) {
  const long MAXENTRYCHARS = 10, MAXCOLS = 1000;
  const long MAXLINELENGTH = MAXCOLS*(MAXENTRYCHARS+1);
  long nrows=0, ncols=0;
  gsl_matrix *table;
  FILE *file, *stream;
  char phrase[MAXLINELENGTH+1], word[MAXENTRYCHARS+1];
  int status;

  // Open file:
  if ((file=fopen(filename.c_str(), "r"))==NULL) error("LoadGSLMatrix: cannot open file "+filename);

  // Count lines and columns:
  fgets(phrase, MAXLINELENGTH, file); nrows++;
  if(strlen(phrase)>MAXLINELENGTH-2) error("LoadGSLMatrix: increase MAXLINELENGTH if possible.");
  stream=fmemopen(phrase, strlen(phrase), "r");
  while (fscanf(stream, "%s", word)!=EOF) ncols++;
  while (fgets(phrase, MAXLINELENGTH, file)!=NULL) nrows++;
  rewind(file);
  
  // Allocate memory:
  std::cout<<"LoadGSLMatrix will allocate "<<nrows<<" lines and "<<ncols<<" columns for file "<<filename<<std::endl;
  table = gsl_matrix_alloc(nrows, ncols);
  if (table==NULL) error("LoadGSLMatrix: gsl_matrix_alloc failed!");

  // Loading values to table:
  status=gsl_matrix_fscanf(file, table);
  if(status==GSL_EFAILED) error("LoadGSLMatrix: gsl_matrix_fscanf failed!");

  fclose(file); fclose(stream);
  return table;
}


/*** Print GSL matrix as a table (with rows and columns) ***/
void PrintGSLMatrix(const gsl_matrix *A, std::ostream *output /*= &std::cout*/) {
  long i, j;

  (*output).setf(std::ios_base::showpoint);
  (*output).precision(14);
  for (i=0; i<(A->size1); i++) {
    for (j=0; j<(A->size2); j++) {
      (*output).width(20); *output << A->data[i*(A->size1)+j] << " ";
    }
    *output << std::endl;
  }       

}
