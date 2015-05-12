#include <healpix_map_fitsio.h> // For FITS files.
#include <levels_facilities.h>  // For something related to Healpix I don't remember.
#include <iomanip>              // For setprecision.
#include "GeneralOutput.hpp"    // I don't know why, maybe to avoid mismatches.
#include "Utilities.hpp"        // For warnings, errros and dynamic allocation.
#include "corrlnfields_aux.hpp" // For n2fz function.

/**********************************************/
/*** Auxiliary functions for General output ***/
/**********************************************/

/*** Function for writing the header of alm's file ***/
std::string SampleHeader(int N1, int N2) {
  std::stringstream ss;
  int i, f, z, nr;
  
  nr=N1*N2;

  ss << "# l, m";
  for (i=0; i<nr; i++) {
    n2fz(i, &f, &z, N1, N2);
    ss << ", f" << f << "z" << z << "Re, f" << f << "z" << z << "Im";
  }
  
  return ss.str(); 
}


/***********************/
/*** Matrices output ***/
/***********************/

// Print a GSL matrix to file
void GeneralOutput(const gsl_matrix *Cov, std::string filename, bool inform) {
  std::ofstream outfile; 

  outfile.open(filename.c_str());
  if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
  else { 
    PrintGSLMatrix(Cov, &outfile); 
    outfile.close();
    if(inform==1) std::cout << ">> GSL matrix written to "+filename<<std::endl;
  }  
}


// Print all GSL matrices in a vector to files:
void GeneralOutput(gsl_matrix **CovByl, const ParameterList & config, std::string keyword, bool inform) {
  std::string filename;
  std::ofstream outfile; 
  int lmin, lmax, l;

  if (config.reads(keyword)!="0") {
    lmin = config.readi("LMIN");
    lmax = config.readi("LMAX");
    for (l=lmin; l<=lmax; l++) {
      filename=config.reads(keyword)+"l"+ZeroPad(l,lmax)+".dat";
      outfile.open(filename.c_str());
      if (!outfile.is_open()) warning("corrlnfields: cannot open file "+filename);
      else { 
	PrintGSLMatrix(CovByl[l], &outfile); 
	outfile.close();
	if(inform==1) std::cout << ">> "+keyword+" ["<<l<<"] written to "+filename<<std::endl;
      }
    }  
  }
}


/********************/
/*** Alm's output ***/
/********************/


// Prints ALL fields alm's to a TEXT file labeled by their FIELD IDs.
void GeneralOutput(Alm<xcomplex <double> > *af, const ParameterList & config, std::string keyword, int N1, int N2, bool inform) {
  std::string filename;
  std::ofstream outfile; 
  int lminout, lmaxout, mmax, l, m, i, Nfields;

  Nfields=N1*N2;

  // If requested, write alm's to file:
  if (config.reads(keyword)!="0") {
    filename = config.reads(keyword);
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("GeneralOutput: cannot open "+filename+" file.");
    outfile << SampleHeader(N1, N2) <<std::endl<<std::endl;
    lminout = config.readi("LRANGE_OUT", 0);
    lmaxout = config.readi("LRANGE_OUT", 1);
    mmax = config.readi("MMAX_OUT");
    if (mmax>lminout) error ("GeneralOutput: current code only allows MMAX_OUT <= LMIN_OUT.");
    // Output all alm's:
    if (mmax<0) {
      for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< af[i](l,m).re<<" "<<std::setprecision(10)<< af[i](l,m).im;
	  outfile<<std::endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
	  for (i=0; i<Nfields; i++) outfile <<" "<<std::setprecision(10)<< af[i](l,m).re<<" "<<std::setprecision(10)<< af[i](l,m).im;
	  outfile<<std::endl;
	}  
    }
    outfile.close();
    if(inform==1) std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}


// Prints one single alm table to a TEXT file using a PREFIX and a FIELD ID.
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, int f, int z, bool inform) {
  std::string filename;
  std::ofstream outfile; 
  char message[100];
  int lminout, lmaxout, mmax, l, m;

  // If requested, write alm's to file:
  if (config.reads(keyword)!="0") {
    sprintf(message, "%sf%dz%d.dat", config.reads(keyword).c_str(), f, z);
    filename.assign(message);
    
    outfile.open(message);
    if (!outfile.is_open()) warning("GeneralOutput: cannot open "+filename+" file.");
    outfile << "# l, m, Re(alm), Im(alm)"<<std::endl<<std::endl;
    lminout = config.readi("LRANGE_OUT", 0);
    lmaxout = config.readi("LRANGE_OUT", 1);
    mmax = config.readi("MMAX_OUT");
    if (mmax>lminout) error ("GeneralOutput: current code only allows MMAX_OUT <= LMIN_OUT.");
    // Output all alm's:
    if (mmax<0) {
      for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
	  outfile<<std::endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
	  outfile<<std::endl;
	}  
    }
    outfile.close();
    if(inform==1) std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}



// Prints one single alm table to a TEXT file.
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, bool inform) {
  std::string filename;
  std::ofstream outfile; 
  int lminout, lmaxout, mmax, l, m;

  // If requested, write alm's to file:
  if (config.reads(keyword)!="0") {
    filename = config.reads(keyword);
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("GeneralOutput: cannot open "+filename+" file.");
    outfile << "# l, m, Re(alm), Im(alm)"<<std::endl<<std::endl;
    lminout = config.readi("LRANGE_OUT", 0);
    lmaxout = config.readi("LRANGE_OUT", 1);
    mmax = config.readi("MMAX_OUT");
    if (mmax>lminout) error ("GeneralOutput: current code only allows MMAX_OUT <= LMIN_OUT.");
    // Output all alm's:
    if (mmax<0) {
      for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=l; m++) {
	  outfile << l <<" "<< m;
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
	  outfile<<std::endl;
	} 
    }
    // Truncate m in alm output:
    else {
     for(l=lminout; l<=lmaxout; l++)
	for(m=0; m<=mmax; m++) {
	  outfile << l <<" "<< m;
	  outfile <<" "<<std::setprecision(10)<< a(l,m).re<<" "<<std::setprecision(10)<< a(l,m).im;
	  outfile<<std::endl;
	}  
    }
    outfile.close();
    if(inform==1) std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}



/********************/
/*** Maps output  ***/
/********************/


// Prints a list of maps to a many FITS files:
void GeneralOutputFITS(Healpix_Map<double> *mapf, const ParameterList & config, std::string keyword, int N1, int N2, bool inform) {
  std::string filename, tempstr;
  char message[100], message2[100], *arg[4];
  char opt1[]="-bar";
  int i, Nfields, f, z;
  
  Nfields=N1*N2;

  if (config.reads(keyword)!="0") {
    // Write to FITS:
    tempstr  = config.reads(keyword);
    for (i=0; i<Nfields; i++) {
      n2fz(i, &f, &z, N1, N2); // Get field and redshift numbers.
      sprintf(message, "%sf%dz%d.fits", tempstr.c_str(), f, z);
      filename.assign(message);
      //sprintf(message2, "rm -f %s", message);
      //system(message2); // Have to delete previous fits files first.
      write_Healpix_map_to_fits("!"+filename,mapf[i],planckType<double>()); // Filename prefixed by ! to overwrite.
      if(inform==1) std::cout << ">> "<<keyword<< "["<<i<<"] written to "<<filename<<std::endl;
      // Write to TGA if requested:
      if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
	sprintf(message2, "%sf%dz%d.tga", tempstr.c_str(), f, z);
	arg[1]=message; arg[2]=message2; arg[3]=opt1;
	map2tga_module(4, (const char **)arg);
	if(inform==1) std::cout << ">> "<<keyword<< "["<<i<<"] written to "<<message2<<std::endl;
	if (config.readi("FITS2TGA")==2) {
	  sprintf(message2, "rm -f %s", message);
	  system(message2);
	  if(inform==1) std::cout << "-- "<<filename<<" file removed."<<std::endl;
	}
      }
    }
  }
}


// Prints a list of maps to a single TEXT file.
void GeneralOutputTEXT(Healpix_Map<double> *mapf, const ParameterList & config, std::string keyword, int N1, int N2, bool inform) {
  std::string filename;
  std::ofstream outfile;
  int i, j, Nfields, field, z, npixels, *nside, n, coordtype;
  pointing coord;

  if (config.reads(keyword)!="0") {
    Nfields   = N1*N2; 
    coordtype = config.readi("ANGULAR_COORD");
    if (coordtype<0 || coordtype>2) warning("GeneralOutputTEXT: unknown ANGULAR_COORD option, keeping Theta & Phi in radians.");

    // Check which maps are allocated and avoid different-size maps.
    j=0;
    nside = vector<int>(0,Nfields-1);
    for (i=0; i<Nfields; i++) {
      nside[i] = mapf[i].Nside(); // Record all maps Nsides.
      if (nside[i]!=0) {          // Look for allocated maps.
	if (j==0) {               // Initialize j to first non-zero Nside.
	  j=nside[i]; n=i;
	}
	else if (nside[i]!=j) error ("GeneralOutput: maps do not have the same number of pixels.");
      }
    }
    npixels=12*j*j;
    
    // Output to file:
    filename=config.reads(keyword);
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("GeneralOutput: cannot open file "+filename);
    else {
      // Set Headers:
      outfile << "# theta, phi";
      for (i=0; i<Nfields; i++) if (nside[i]!=0) { 
	  n2fz(i, &field, &z, N1, N2);
	  outfile << ", f"<<field<<"z"<<z;
	}
      outfile << std::endl;
      // LOOP over allocated fields:
      for (j=0; j<npixels; j++) {
	coord = mapf[n].pix2ang(j);
	if      (coordtype==1) {coord.theta=rad2deg  (coord.theta); coord.phi=rad2deg(coord.phi);}
	else if (coordtype==2) {coord.theta=theta2dec(coord.theta); coord.phi=phi2ra (coord.phi);}
	outfile << coord.theta <<" "<< coord.phi;
	for (i=0; i<Nfields; i++) if (nside[i]!=0) outfile <<" "<< mapf[i][j];
	outfile << std::endl;
      }
    }
    outfile.close();
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    free_vector(nside, 0,Nfields-1);
  } 
}


// Prints array of Healpix maps to either a single TEXT (if fits=0 or unespecified) or many FITS (if fits=1) files.
void GeneralOutput(Healpix_Map<double> *mapf, const ParameterList & config, std::string keyword, int N1, int N2, bool fits, bool inform) {
  if (fits==1) GeneralOutputFITS(mapf, config, keyword, N1, N2, inform);
  else GeneralOutputTEXT(mapf, config, keyword, N1, N2, inform);
}


// Prints two lists of maps to a single TEXT file.
void GeneralOutput(Healpix_Map<double> *gamma1, Healpix_Map<double> *gamma2, 
		   const ParameterList & config, std::string keyword, int N1, int N2, bool inform) {
  std::string filename;
  std::ofstream outfile;
  int i, j, Nfields, field, z, npixels, *nside, n;
  pointing coord;

  if (config.reads(keyword)!="0") {
    Nfields=N1*N2; 
    
    // Check which maps are allocated and avoid different-size maps.
    j=0;
    nside = vector<int>(0,Nfields-1);
    for (i=0; i<Nfields; i++) {
      nside[i]  = gamma1[i].Nside(); // Record all maps Nsides.
      if (nside[i]!=gamma2[i].Nside()) error ("GeneralOutput: expecting pair of maps with same number of pixels.");
      if (nside[i]!=0) {             // Look for allocated maps.
	if (j==0) {                  // Initialize j to first non-zero Nside.
	  j=nside[i]; n=i;
	}
	else if (nside[i]!=j) error ("GeneralOutput: maps do not have the same number of pixels.");
      }
    }
    npixels=12*j*j;
    
    // Output to file:
    filename=config.reads(keyword);
    outfile.open(filename.c_str());
    if (!outfile.is_open()) warning("GeneralOutput: cannot open file "+filename);
    else {
      // Set Headers:
      outfile << "# theta, phi";
      for (i=0; i<Nfields; i++) if (nside[i]!=0) { 
	  n2fz(i, &field, &z, N1, N2);
	  outfile <<", f"<<field<<"z"<<z<<"[1]"<<", f"<<field<<"z"<<z<<"[2]";
	}
      outfile << std::endl;
      // LOOP over allocated fields:
      for (j=0; j<npixels; j++) {
	coord = gamma1[n].pix2ang(j);
	outfile << coord.theta <<" "<< coord.phi;
	for (i=0; i<Nfields; i++) if (nside[i]!=0) outfile <<" "<< gamma1[i][j]<<" "<< gamma2[i][j];
	outfile << std::endl;
      }
    }
    outfile.close();
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    free_vector(nside, 0,Nfields-1);
  } 
}


// Prints one single map to FITS file based on a PREFIX and a FIELD ID:
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword, int *fnz, bool inform) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    sprintf(message1, "%sf%dz%d.fits", config.reads(keyword).c_str(), fnz[0], fnz[1]);
    filename.assign(message1);

    // Write to FITS:
    //sprintf(message1, "rm -f %s", filename.c_str());
    //system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits("!"+filename, map, planckType<double>()); // Filename prefixed by ! to overwrite.
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      if(inform==1) std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	if(inform==1) std::cout << "-- "<<filename<<" file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single combination of kappa, gamma1 and gamma2 maps to FITS file based on a PREFIX and a FIELD ID:
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, int f, int z, bool inform) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    sprintf(message1, "%sf%dz%d.fits", config.reads(keyword).c_str(), f, z);
    filename.assign(message1);

    // Write to FITS:
    //sprintf(message1, "rm -f %s", filename.c_str());
    //system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits("!"+filename, kmap, g1map, g2map, planckType<double>());
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      if(inform==1) std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	if(inform==1) std::cout << "-- "<<filename<<" file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single combination of kappa, gamma1 and gamma2 maps to FITS file.
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, bool inform) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    // Write to FITS:
    filename=config.reads(keyword); 
    //sprintf(message1, "rm -f %s", filename.c_str());
    //system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits("!"+filename, kmap, g1map, g2map, planckType<double>()); // Filename prefixed by ! to overwrite.
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      if(inform==1) std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	if(inform==1) std::cout << "-- "<<keyword<<" FITS file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single map to FITS and/or TGA file.
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword, bool inform) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    // Write to FITS:
    filename=config.reads(keyword); 
    //sprintf(message1, "rm -f %s", filename.c_str());
    //system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits("!"+filename, map, planckType<double>()); // Filename prefixed by ! to overwrite.
    if(inform==1) std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      if(inform==1) std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	if(inform==1) std::cout << "-- "<<keyword<<" FITS file removed."<<std::endl;
      }
    }
  } 
}
