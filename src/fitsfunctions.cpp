#include <cstring>
#include <fitsio.h>
#include "../src/definitions.hpp"
#include "../src/Utilities.hpp"
#include "../src/flask_aux.hpp"

// The definition below will be changed by make and the result will be written to a new file:
#define HEALPIX_DATA "/path/to/healpix/dir/data"

int WriteCatalog2Fits(std::string filename, CAT_PRECISION **table, long Nentries, const ParameterList & config) {
  const int COLNAMELENGTH=20;
  fitsfile *fpointer;
  std::stringstream ss;
  std::string header, word;
  char **columnNames, **columnTypes, **columnUnits;
  int status=0, Ncols, i, datatype, *integer, fits_precision;
  long j;
  char TableName[]  ="GALAXY_CAT"; 

  // Allocate auxiliary integer vector:
  integer = vector<int>(0, Nentries-1);

  // Get number of columns:
  header      = config.reads("CATALOG_COLS");
  Ncols       = CountWords(header);
  // Reaname 'theta' and 'phi' to 'dec' and 'ra' if change of coords. was requested:
  if (config.readi("ANGULAR_COORD")==2) {
    StrReplace(header, "theta", "dec");
    StrReplace(header, "phi", "ra");
  }
  // Get column names:
  columnNames = matrix<char>(0,Ncols, 0,COLNAMELENGTH);
  ss << header;
  for(i=0; i<Ncols; i++) { ss >> word; strcpy(columnNames[i], word.c_str()); } 
  
  // Set column formats and units:
  columnTypes = matrix<char>(0,Ncols, 0,6);
  columnUnits = matrix<char>(0,Ncols, 0,COLNAMELENGTH);
  for(i=0; i<Ncols; i++) { 
    // theta phi z galtype kappa gamma1 gamma2 ellip1 ellip2 pixel maskbit
    if      (strcmp(columnNames[i],"theta"  )==0) {sprintf(columnTypes[i],"%s", "F9.5" ); sprintf(columnUnits[i],"%s", "Radians");}
    else if (strcmp(columnNames[i],"phi"    )==0) {sprintf(columnTypes[i],"%s", "F9.5" ); sprintf(columnUnits[i],"%s", "Radians");}
    else if (strcmp(columnNames[i],"ra"     )==0) {sprintf(columnTypes[i],"%s", "F9.5" ); sprintf(columnUnits[i],"%s", "Degrees");}
    else if (strcmp(columnNames[i],"dec"    )==0) {sprintf(columnTypes[i],"%s", "F9.5" ); sprintf(columnUnits[i],"%s", "Degrees");}
    else if (strcmp(columnNames[i],"z"      )==0) {sprintf(columnTypes[i],"%s", "F8.5" ); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"galtype")==0) {sprintf(columnTypes[i],"%s", "I3"   ); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"kappa"  )==0) {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"gamma1" )==0) {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"gamma2" )==0) {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"ellip1" )==0) {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"ellip2" )==0) {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"pixel"  )==0) {sprintf(columnTypes[i],"%s", "I8"   ); sprintf(columnUnits[i],"%s", "\0");}
    else if (strcmp(columnNames[i],"maskbit")==0) {sprintf(columnTypes[i],"%s", "I2"   ); sprintf(columnUnits[i],"%s", "\0");}
    else                                          {sprintf(columnTypes[i],"%s", "E12.5"); sprintf(columnUnits[i],"%s", "unknown"); 
      warning("WriteCatalog2Fits: unknown catalog column "+word.assign(columnNames[i]));}  
  }

  // Create (or overwrite) FITS file with ASCII table:
  fits_create_file(&fpointer, ("!"+filename).c_str(), &status);
  fits_report_error(stderr, status);
  fits_create_tbl(fpointer, ASCII_TBL, Nentries, Ncols, columnNames, columnTypes, columnUnits, TableName, &status);
  fits_report_error(stderr, status);
  

  // Write columns to FITS file ASCII table:
  for(i=0; i<Ncols; i++) {
    // Write double variable:
    if (columnTypes[i][0] == 'F' || columnTypes[i][0] == 'E' || columnTypes[i][0] == 'D') {
      fits_write_col(fpointer, FIT_PRECISION, i+1, 1, 0, Nentries, table[i], &status);
    }
    // Write int variable:
    else if (columnTypes[i][0] == 'I') {
      for(j=0; j<Nentries; j++) integer[j]=(int)table[i][j];
      fits_write_col(fpointer, TINT, i+1, 1, 0, Nentries, integer, &status);
    }
    else error("WriteCatalog2Fits: "+word.assign(columnNames[i])+" has uknown FITS format.");
    // Verify if everything is OK:
    fits_report_error(stderr, status);
  }
  
  // Do not need columns info anymore:
  free_matrix(columnNames, 0,Ncols, 0,COLNAMELENGTH);
  free_matrix(columnTypes, 0,Ncols, 0,6);
  free_matrix(columnUnits, 0,Ncols, 0,COLNAMELENGTH);
  // Free auxiliary integer vector:
  free_vector(integer, 0,Nentries-1);

  // Close FITS file and exit:
  fits_close_file(fpointer, &status);
  fits_report_error(stderr, status);
  return status;
}


// Reads a Healpix FITS file containing the map2alm weights into a double array:
int ReadHealpixData(int col, const ParameterList & config, double *weights, int datatype) {
  char message[200];
  std::string filename;
  fitsfile *fpointer;
  int status=0, anynul=0;
  long i, firstrow=1, firstelem=1, nelements, nside;
  double *nulval;
  
  nside    = config.readi("NSIDE");
  filename.assign(HEALPIX_DATA);
  if (filename.at(filename.length()-1)!='/') filename = filename+"/";
  // Find out what kind of data will be loaded:
  switch (datatype) {
  case 1:
    // Load map2alm weights:
    filename = filename+"weight_ring_n"+ZeroPad(nside, 10000)+".fits";
    nelements=2*nside;
    break;
  case 2:
    // Load pixel window functions:
    filename = filename+"pixel_window_n"+ZeroPad(nside, 1000)+".fits";
    nelements=4*nside+1;
    break;
  default:
    error("ReadHealpixData: unknown Healpix data type.");
  }

  // Open file:
  fits_open_table(&fpointer, filename.c_str(), READONLY, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: could not open table in FITS file, ERR=%d", status);
    warning(message);
  }
  
  // Prepare to, read and check for errors:
  nulval = vector<double>(0, nelements-1);
  for(i=0; i<nelements; i++) nulval[i]=666.0;
  fits_read_col(fpointer, TDOUBLE, col, firstrow, firstelem, nelements, nulval, weights, &anynul, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: problem reading column in FITS file table, ERR=%d", status);
    warning(message);
  }
  if(anynul!=0) {
    warning("ReadHealpixData: found NULL values in FITS table");
    printf("They are:\n");
    for (i=0; i<nelements; i++) printf("%g ",nulval[i]);
    printf("\n");
  }
  free_vector(nulval, 0, nelements-1);

  // Close file and exit:
  fits_close_file(fpointer, &status);
  if (status!=0) {
    sprintf(message, "ReadHealpixData: could not close FITS file, ERR=%d", status);
    warning(message);
  }
  return status;
}
