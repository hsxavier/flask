#include <healpix_map_fitsio.h> // For GeneralOutput function.
#include <levels_facilities.h>  // For GeneralOutput function.
#include <iomanip>              // For GeneralOutput function.
#include "GeneralOutput.hpp"

// Prints one single combination of kappa, gamma1 and gamma2 maps to FITS file based on a PREFIX and a FIELD ID:
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword, int *fnz) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    sprintf(message1, "%sf%dz%d.fits", config.reads(keyword).c_str(), fnz[0], fnz[1]);
    filename.assign(message1);

    // Write to FITS:
    sprintf(message1, "rm -f %s", filename.c_str());
    system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits(filename, kmap, g1map, g2map, planckType<double>());
    std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	std::cout << "-- "<<filename<<" file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single combination of kappa, gamma1 and gamma2 maps to FITS file.
void GeneralOutput(const Healpix_Map<double> & kmap, const Healpix_Map<double> & g1map, 
		   const Healpix_Map<double> & g2map, const ParameterList & config, std::string keyword) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    // Write to FITS:
    filename=config.reads(keyword);
    sprintf(message1, "rm -f %s", filename.c_str());
    system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits(filename, kmap, g1map, g2map, planckType<double>());
    std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	std::cout << "-- "<<keyword<<" FITS file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single map to FITS and/or TGA file.
void GeneralOutput(const Healpix_Map<double> & map, const ParameterList & config, std::string keyword) {
  std::string filename, tgafile;
  char *arg[4];
  char message1[100], message2[100];
  char opt1[]="-bar";
  if (config.reads(keyword)!="0") {
    // Write to FITS:
    filename=config.reads(keyword);
    sprintf(message1, "rm -f %s", filename.c_str());
    system(message1); // Have to delete previous fits files first.
    write_Healpix_map_to_fits(filename, map, planckType<double>());
    std::cout << ">> "<<keyword<<" written to "<<filename<<std::endl;
    // Write to TGA if requested:
    if (config.readi("FITS2TGA")==1 || config.readi("FITS2TGA")==2) {
      tgafile = filename;
      tgafile.replace(tgafile.find(".fits"),5,".tga");
      sprintf(message1, "%s", filename.c_str());
      sprintf(message2, "%s", tgafile.c_str());
      arg[1]=message1; arg[2]=message2; arg[3]=opt1;
      map2tga_module(4, (const char **)arg);
      std::cout << ">> "<<keyword<<" written to "<<tgafile<<std::endl;
      if (config.readi("FITS2TGA")==2) {
	sprintf(message2, "rm -f %s", message1);
	system(message2);
	std::cout << "-- "<<keyword<<" FITS file removed."<<std::endl;
      }
    }
  } 
}


// Prints one single alm table to a TEXT file using a PREFIX and a FIELD ID.
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword, int *fnz) {
  std::string filename;
  std::ofstream outfile; 
  char message[100];
  int lminout, lmaxout, mmax, l, m;

  // If requested, write alm's to file:
  if (config.reads(keyword)!="0") {
    sprintf(message, "%sf%dz%d.dat", config.reads(keyword).c_str(), fnz[0], fnz[1]);
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
    std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}



// Prints one single alm table to a TEXT file.
void GeneralOutput(const Alm<xcomplex <double> > & a, const ParameterList & config, std::string keyword) {
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
    std::cout << ">> "+keyword+" written to "+filename<<std::endl;
  }
}
