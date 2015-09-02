// These are functions created specifically for the corrlnfields program.
#include "corrlnfields_aux.hpp"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <xcomplex.h>
#include "fitsfunctions.hpp"
#include "GeneralOutput.hpp"
#include <alm_healpix_tools.h>
#include "lognormal.hpp" // For gmu, gsigma, etc. in PrintMapsStats function.


bool ComputeShearQ(const ParameterList & config) {
  std::string CatalogHeader, ExitAt;
  
  // Only check if output is required if code does not stop before them:
  ExitAt = config.reads("EXIT_AT");
  if (ExitAt!="FLIST_OUT"       && ExitAt!="XIOUT_PREFIX"    && ExitAt!="GXIOUT_PREFIX"   && 
      ExitAt!="GCLOUT_PREFIX"   && ExitAt!="COVL_PREFIX"     && ExitAt!="REG_COVL_PREFIX" && 
      ExitAt!="REG_CL_PREFIX"   && ExitAt!="CHOLESKY_PREFIX" && ExitAt!="AUXALM_OUT"      && 
      ExitAt!="AUXMAP_OUT"      && ExitAt!="MAP_OUT"         && ExitAt!="MAPFITS_PREFIX"  && 
      ExitAt!="DENS2KAPPA_STAT" && ExitAt!="RECOVALM_OUT"    && ExitAt!="RECOVCLS_OUT"    ) {

    // Return 1 if any shear output is required:
    if (config.reads("SHEAR_ALM_PREFIX")  !="0") return 1;
    if (config.reads("SHEAR_FITS_PREFIX") !="0") return 1;
    if (config.reads("SHEAR_MAP_OUT")     !="0") return 1;
    if (config.reads("CATALOG_OUT")       !="0") {
      CatalogHeader = config.reads("CATALOG_COLS");
      if (GetSubstrPos("gamma1" , CatalogHeader) != -1) return 1;
      if (GetSubstrPos("gamma2" , CatalogHeader) != -1) return 1;
      if (GetSubstrPos("ellip1" , CatalogHeader) != -1) return 1;
      if (GetSubstrPos("ellip2" , CatalogHeader) != -1) return 1;
    }
  
  } // End of IF stop code before shear output.

  // If no shear related output is requested, return 0:
  return 0;
}


double MapMean(const Healpix_Map<MAP_PRECISION> & map) {
  return map.average();
}

double MapVariance(const Healpix_Map<MAP_PRECISION> & map, double mean) {
  int i, npixels;
  double var=0.0, aux;
  
  npixels = map.Npix();
#pragma omp parallel for reduction(+:var) private(aux)
  for (i=0; i<npixels; i++) {
    aux  = map[i] - mean;
    var += aux*aux;
  }
  return var/npixels;
}

double MapSkewness(const Healpix_Map<MAP_PRECISION> & map, double mean, double variance) {
  int i, npixels;
  double skew=0.0, aux;
  
  npixels = map.Npix();
#pragma omp parallel for reduction(+:skew) private(aux)
  for (i=0; i<npixels; i++) {
    aux  = map[i] - mean;
    skew += pow(aux, 3);
  }
  return skew/npixels/pow(variance, 1.5);
}


// Prints a table with mean, std. dev., skewness, mu (gaussian mean), sigma (gaussian dev) and shift for all maps.
void PrintMapsStats(Healpix_Map<MAP_PRECISION> *mapf, const FZdatabase & fieldlist, simtype dist, std::ostream *output) {
  const int ColWidth = 12;
  int Nfields, i, f, z;
  char name[11];
  double mean, var, skew, shift;
  
  Nfields = fieldlist.Nfields();
  (*output).width(ColWidth);
  
  // Header:
  (*output) <<std::left<<"# FieldID"<<std::right; (*output).width(ColWidth);
  (*output) << "Mean";           (*output).width(ColWidth);
  (*output) << "Std.Dev.";       (*output).width(ColWidth);
  (*output) << "Skewness";       (*output).width(ColWidth);
  if (dist==lognormal) {
    (*output) << "gMU";          (*output).width(ColWidth);
    (*output) << "gSIGMA";       (*output).width(ColWidth);
    (*output) << "Shift";        (*output).width(ColWidth);
  }
  (*output) << std::endl;        (*output).width(ColWidth);

  // LOOP over allocated maps: 
  for (i=0; i<Nfields; i++) if (mapf[i].Nside()>0) {
      // Compute statistics:
      mean  = MapMean(mapf[i]);
      var   = MapVariance(mapf[i], mean);
      skew  = MapSkewness(mapf[i], mean, var);
      shift = Moments2Shift(mean, var, skew); 
      // Output results:
      fieldlist.Index2Name(i, &f, &z);
      sprintf(name, "f%dz%d   ",f,z);
      (*output) << std::left<<name<<std::right; (*output).width(ColWidth);
      (*output) << mean;                        (*output).width(ColWidth);
      (*output) << sqrt(var);                   (*output).width(ColWidth);
      (*output) << skew;                        (*output).width(ColWidth);
      if (dist==lognormal) {
	(*output) << gmu(mean, var, shift);     (*output).width(ColWidth);
	(*output) << gsigma(mean, var, shift);  (*output).width(ColWidth);
	(*output) << shift;                     (*output).width(ColWidth);
      }
      (*output) << std::endl;                   (*output).width(ColWidth);
    }
}

// If alm or Cls output is requested, compute them and output them:
void RecoverAlmCls(Healpix_Map<MAP_PRECISION> *mapf, const FZdatabase & fieldlist, 
		   std::string almKey, std::string clsKey, const ParameterList & config) {
  int i, j, k, l, m, nside, lmin, lmax, lminout, lmaxout, mmax, NCls, Nfields;
  Alm<xcomplex <ALM_PRECISION> > *bflm;
  double **recovCl;
  bool *yesCl;

  if (config.reads(almKey)!="0" || config.reads(clsKey)!="0") {

    Nfields = fieldlist.Nfields();
    nside   = config.readi("NSIDE");
    lmin    = config.readi("LMIN");
    lmax    = config.readi("LMAX");  
    

    // Allocate and clear variables to receive new alm's:
    bflm  = vector<Alm<xcomplex <ALM_PRECISION> > >(0,Nfields-1);
    for (i=0; i<Nfields; i++) if(mapf[i].Nside()!=0) {
	bflm[i].Set(lmax,lmax);
	for(l=0; l<=lmax; l++) for (m=0; m<=l; m++) bflm[i](l,m).Set(0,0);
      }    
    // Load map ring weights:
    arr<double> weight(2*nside);
    PrepRingWeights(1, weight, config);
    // Compute alm's from map:
    Announce("Recovering alm's from map... ");
    for(i=0; i<Nfields; i++) if(mapf[i].Nside()!=0) map2alm(mapf[i],bflm[i],weight);
    Announce();
    weight.dealloc();
    // Output to file:
    GeneralOutput(bflm, config, almKey, fieldlist);

    // Compute Cl's if requested:
    if (config.reads(clsKey)!="0") {
      Announce("Recovering Cl's from maps... ");
      lminout = config.readi("LRANGE_OUT", 0);
      lmaxout = config.readi("LRANGE_OUT", 1);
      if (lmaxout > lmax) { lmaxout=lmax; warning("RecoverAlmCls: LRANGE_OUT beyond LMAX, will use LMAX instead"); }
      if (lminout < lmin) { lminout=lmin; warning("RecoverAlmCls: LRANGE_OUT beyond LMIN, will use LMIN instead"); }
      mmax    = config.readi("MMAX_OUT");
      if (mmax>lminout) error("RecoverAlmCls: current code only allows MMAX_OUT <= LMIN_OUT.");
      NCls  = Nfields*(Nfields+1)/2;
      recovCl = matrix<double>(0, NCls-1, lminout, lmaxout);
      yesCl   = vector<bool>(0, NCls-1);
#pragma omp parallel for schedule(dynamic) private(l, m, i, j)
      for (k=0; k<NCls; k++) {
	l = (int)((sqrt(8.0*(NCls-1-k)+1.0)-1.0)/2.0);
	m = NCls-1-k-(l*(l+1))/2;
	i = Nfields-1-l;
	j = Nfields-1-m;
	if(bflm[i].Lmax()!=0 && bflm[j].Lmax()!=0) {
	  yesCl[k] = 1;
	  for(l=lminout; l<=lmaxout; l++) {
	    recovCl[k][l] = 0;
	    if (mmax<0) for (m=0; m<=l; m++)    recovCl[k][l] += (bflm[i](l,m)*(bflm[j](l,m).conj())).real();
	    else        for (m=0; m<=mmax; m++) recovCl[k][l] += (bflm[i](l,m)*(bflm[j](l,m).conj())).real();
	    recovCl[k][l] = recovCl[k][l]/((double)(l+1));
	  }
	}
	else yesCl[k] = 0;
      } // End of Cl computing.
      Announce();
      GeneralOutput(recovCl, yesCl, fieldlist, config, clsKey);
      free_matrix(recovCl, 0, NCls-1, lminout, lmaxout);
      free_vector(yesCl, 0, NCls-1);
    }
    
    // Free memory:
    free_vector(bflm, 0, Nfields-1);
  }
}


// Prepare weights used by map2alm. weight array must be allocated already:
void PrepRingWeights(int col, arr<double> & weight, const ParameterList & config) {
  double *tempweight;
  int nside, status;
  long i;

  nside = config.readi("NSIDE");

  if (config.readi("USE_HEALPIX_WGTS")==1) {
    Announce("   Loading Healpix map weights... ");
    tempweight = vector<double>(0, 2*nside-1);
    status = ReadHealpixWeights(col, nside, config, tempweight);
    if (status==0) for (i=0; i<2*nside; i++) weight[i]=1.0+tempweight[i];
    else { 
      warning("PrepRingWeights: could not load Healpix weights, using 1.0 instead.");
      weight.fill(1);
    }
    free_vector(tempweight, 0, 2*nside-1);
    Announce();
  }
  else weight.fill(1);  
}


// Transform radians to degrees:
double rad2deg(double rad) {
  const double OneEightyOverPi = 180.0/M_PI;
  return OneEightyOverPi * rad;
}

// Transform Theta (in radians) into DEC (in degrees):
double theta2dec(double theta) {
  const double PiOverTwo = M_PI/2.0;
  return rad2deg(PiOverTwo - theta);
}

// Transform Phi (in radians) into RA (in degrees):
double phi2ra(double phi) {
  return rad2deg(phi);
}


// Get file format from filename extension: .dat, .fits, .fits.gz
int FileFormat(std::string filename) {
  using namespace definitions;
  int pos;
  
  // Possible file extensions:
  pos = filename.rfind(".fits.gz");
  if (pos>0 && pos+8==filename.size()) return fits_format;
  pos = filename.rfind(".fits");
  if (pos>0 && pos+5==filename.size()) return fits_format;
  pos = filename.rfind(".dat");
  if (pos>0 && pos+4==filename.size()) return ascii_format;
  
  // Unidentified format:
  return unknown_format;
}


// Count words in phrase:
int CountWords(const std::string header) {
  std::string entry;
  std::stringstream ss(header);
  int column=0;

  while (ss >> entry) column++;
  return column;
}


// Return position of word in phrase, return -1 if word is not found:
int GetSubstrPos(const std::string field, const std::string header) {
  std::string entry;
  std::stringstream ss(header);
  int column=0;

  while (ss >> entry && entry!=field) column++;
  if (entry!=field) return -1;
  return column;
}


// Unless column=-1, write value to catalog[column][row] and update catSet:
// Note that the catalog is transposed to ease FITS outputting.
void CatalogFill(CAT_PRECISION **catalog, long row, int column, double value, char **catSet) {  
  if(column==-1) return;
  
  // Write to catalog:            v Count number of updates in cell (for bookkeeping).
  catalog[column][row] = value;   catSet[column][row]++;
}


// Get the shear E-mode harmonic coefficients from the convergence harmonic coefficients:
// (can be done in place)
void Kappa2ShearEmode(Alm<xcomplex <ALM_PRECISION> > &Elm, Alm<xcomplex <ALM_PRECISION> > &Klm) {
  int l, m, lmax;
  double coeff;

  if (Elm.Lmax()!=Klm.Lmax()) error("Kappa2ShearEmode: Elm and klm must have the same lmax.");
  if (Elm.Mmax()!=Klm.Mmax()) error("Kappa2ShearEmode: Elm and klm must have the same mmax.");  
  lmax = Klm.Lmax();
  
  // E mode monopole and dipole are zero:
  for(l=0; l<2; l++) for (m=0; m<=l; m++) Elm(l,m).Set(0,0);

  // Use Wayne Hu (PRD 62:043007, 2000) to get Elm from klm:
  for(l=2; l<=lmax; l++) { 
    coeff = sqrt( ((double)((l+2)*(l-1))) / ((double)(l*(l+1))) );
    for (m=0; m<=l; m++) Elm(l,m).Set(coeff*Klm(l,m).re,coeff*Klm(l,m).im);
  }
}


// Generates galaxy ellipticity from shear and convergence, including random source ellipticity:
void GenEllip(gsl_rng *rnd, double sigma, double kappa, double gamma1, double gamma2, double *eps1, double *eps2) {
  xcomplex<double> g, epsSrc, eps, k, one, gamma;

  // Set complex numbers:
  gamma.re = gamma1; gamma.im = gamma2;
  k.re     = kappa;  k.im     = 0.0;
  one.re   = 1.0;    one.im   = 0.0;
  // Compute reduced shear:
  g = gamma/(one-k);
  // Generate source intrinsic ellipticity:
  if (sigma>0.0) {
    epsSrc.re = gsl_ran_gaussian(rnd, sigma);
    epsSrc.im = gsl_ran_gaussian(rnd, sigma);
   }
  else {
    epsSrc.re = 0.0;
    epsSrc.im = 0.0;
  }
  // Compute ellipticity of the image:
  if (g.norm() <= 1.0) {
    eps = (epsSrc+g) / (one + g.conj()*epsSrc);
  }
  else {
    eps = (one + g*epsSrc.conj())/(epsSrc.conj()+g.conj());
  }
  //Return ellipticity of the image:
  (*eps1) = eps.re;
  (*eps2) = eps.im;
} 


// Uniformly randomly selects an angular position inside a pixel.
pointing RandAngInPix(gsl_rng *r, const Healpix_Map<MAP_PRECISION> & map, int pixel) {
  const double twopi=6.283185307179586;
  std::vector<vec3> corner;
  double thetamin, thetamax, phimin, phimax;
  pointing ang;
  // Find pixel limits for random sampling the angles:
  map.boundaries(pixel, 1, corner);
  thetamin = xyz2ang(corner[0]).theta; // N corner.
  thetamax = xyz2ang(corner[2]).theta; // S corner.
  phimin   = xyz2ang(corner[1]).phi;   // W corner.
  phimax   = xyz2ang(corner[3]).phi;   // E corner.
  if (phimin>phimax) phimin = phimin-twopi;
  // Randomly pick angle inside the pixel. 
  do {ang = randang(r, thetamin, thetamax, phimin, phimax);} 
  while (map.ang2pix(ang) != pixel);
  // Return pointing:
  return ang;
}

// Randomly picks an angular position (uniformly) inside angular boundaries. Needed for RandAngInPix. 
pointing randang(gsl_rng *r, double thetamin, double thetamax, double phimin, double phimax) {
  const double twopi=6.283185307179586;
  double xmin, xmax;
  pointing ang;

  xmin      = (1.0+cos(thetamax))/2.0;
  xmax      = (1.0+cos(thetamin))/2.0;
  ang.phi   = gsl_rng_uniform(r)*(phimax-phimin)+phimin;
  ang.theta = acos(2*(gsl_rng_uniform(r)*(xmax-xmin)+xmin)-1);
  return ang;
}

// Transforms normalized {x,y,z} cartesian coordinates to unit spherical {theta, phi}:
pointing xyz2ang(const vec3 & cartesian) {
  const double pi=3.141592653589793;
  double tempphi;
  pointing ang;
  ang.theta = acos(cartesian.z);
  tempphi   = atan(cartesian.y/cartesian.x);
  if   (cartesian.x<0)    ang.phi = tempphi + pi;   // Quadrants 2 & 3.
  else {
    if (cartesian.y>0)    ang.phi = tempphi;        // Quadrant  1.
    else                  ang.phi = tempphi + pi*2; // Quadrant  4.
  }
  return ang;
}

// Returns the coefficients (x,y,z) of a vector described in a basis that is rotated by 'ang' from the original basis. 
vec3 VecInRotBasis(const pointing & ang, const vec3 & orig) {
  vec3 nuovo;
  nuovo.x = cos(ang.theta)*cos(ang.phi)*orig.x + cos(ang.theta)*sin(ang.phi)*orig.y - sin(ang.theta)*orig.z;
  nuovo.y =      -1.0*sin(ang.phi)     *orig.x +          cos(ang.phi)      *orig.y;
  nuovo.z = sin(ang.theta)*cos(ang.phi)*orig.x + sin(ang.theta)*sin(ang.phi)*orig.y + cos(ang.theta)*orig.z;
  return nuovo;
}

// Uniformly samples the redshift of a bin. THIS IS WRONG, OR AN APPROXIMATION FOR VERY FINE BINS.
double RandRedshift0(gsl_rng *r, double zmin, double zmax) {
  return gsl_rng_uniform(r)*(zmax-zmin)+zmin;
}

// Wrapper so ProjDensityIntegrand can be used by GSL minimizer:
double gsl_mProjDensityIntegrand(double z, void *p) {
  Cosmology *params;
  params = (Cosmology*)p;
  return -1.0*ProjDensityIntegrand(z, params);
}

// Selection of random redshift inside a bin. It is wrong since the density should 
// take into account the selection function, which currently it does not. 
// So it is easier to sample uniformly and make the bins small.
// This code was left here in case we need to implement a GSL minimizer later on
// or something like that.
double ran_redshift(gsl_rng *r, double zmin, double zmax, Cosmology *p) {
  const int MAXITER=100;
  const  double ymin=0, zprecision=0.0005;
  static bool   init=0;
  static double YMAX;
  double ymax, zguess;
  int status, iter=0;
  
  // Find global minimum to initialize the function:
  if (init==0) {
    // Initialize the GSL minimizer:
    gsl_function gslF;
    gslF.function = &gsl_mProjDensityIntegrand;
    gslF.params   = p;
    const gsl_min_fminimizer_type * mintype = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *minscheme = gsl_min_fminimizer_alloc(mintype);
    zguess = (zmax-zmin)/2.0;
    status = gsl_min_fminimizer_set(minscheme, &gslF, zguess, zmin, zmax);
    if (status!=0) error("ran_redshift: problems with gsl_min_fminimizer_set.");
    // Run the minimizer:
    do {
      iter++;
      status = gsl_min_fminimizer_iterate   (minscheme);
      zguess = gsl_min_fminimizer_x_minimum (minscheme);
      zmin   = gsl_min_fminimizer_x_lower   (minscheme);
      zmax   = gsl_min_fminimizer_x_upper   (minscheme);
      status = gsl_min_test_interval(zmin, zmax, zprecision, 0.0);
      if (status == GSL_SUCCESS) printf ("Converged:\n");
      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, zmin, zmax, zguess, zguess - 2.5, zmax - zmin);
    } while (status == GSL_CONTINUE && iter < MAXITER);
    gsl_min_fminimizer_free(minscheme);
    init=1;
    return zguess;
  }

  
  gsl_rng_uniform(r);
}

/*** Multiply Lower-triangular matrix L to complex vector gaus0 and return gaus1 ***/
void CorrGauss(double **gaus1, gsl_matrix *L, double **gaus0) {
  long i, j;

  for(i=0; i<L->size1; i++) {
    gaus1[i][0]=0; // Re    
    gaus1[i][1]=0; // Im
    for(j=0; j<=i; j++) { // L matrix stored as vector in row-major order.
      gaus1[i][0] += L->data[i*L->size1+j] * gaus0[j][0];
      gaus1[i][1] += L->data[i*L->size1+j] * gaus0[j][1];
    }
  }
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


/*** Assign a matrix column n to a variable 'a' identified by a1 and a2  ***/
void fz2n (int a1, int a2, int *n, int N1, int N2) {
  if (a2>N2 || a1>N1 || a1<1 || a2<1) warning("fz2n: unexpected input values.");
  *n = (a1-1)*N2+a2-1; 
}


/*** The inverse of ij2fzfz above ***/
void n2fz (int n, int *a1, int *a2, int N1, int N2) {
  if (n<0 || n>=N1*N2) warning("n2fz: unexpected input values.");
  *a2 = n%N2+1;
  *a1 = n/N2+1;
}


/*** Assign a matrix row i to a variable 'a' identified by a1 and a2 ***/
/*** Assign a matrix column j to a variable 'b' identified by b1 and b2  ***/
void fzfz2ij (int a1, int a2, int b1, int b2, int *i, int *j, int N1, int N2) {
  if (a2>N2 || b2>N2 || a1>N1 || b1>N1 || a1<1 || a2<1 || b1<1 || b2<1) warning("fzfz2ij: unexpected input values.");
  fz2n(a1, a2, i, N1, N2);
  fz2n(b1, b2, j, N1, N2);
}

/*** The inverse of ij2fzfz above ***/
void ij2fzfz (int i, int j, int *a1, int *a2, int *b1, int *b2, int N1, int N2) {
  if (i<0 || j<0 || i>=N1*N2 || j>=N1*N2) warning("ij2fzfz: unexpected input values.");
  n2fz(i, a1, a2, N1, N2);
  n2fz(j, b1, b2, N1, N2);
}


/*** Function for testing the assignments above ***/
void test_fzij (int N1, int N2) {
  int a1, a2, b1, b2, i, j, newa1, newa2, newb1, newb2;
  bool **IsSet;

  IsSet = matrix<bool>(0,N1*N2-1,0,N1*N2-1);
  for(i=0; i<N1*N2; i++) for(j=0; j<N1*N2; j++) IsSet[i][j]=0;
  
  for (a1=1; a1<=N1; a1++)
    for (a2=1; a2<=N2; a2++)
      for (b1=1; b1<=N1; b1++)
	for (b2=1; b2<=N2; b2++) {
	  fzfz2ij(a1, a2, b1, b2, &i, &j, N1, N2); 
	  if (IsSet[i][j]==1) error("test_fzij: tried to set [i,j] already set.");
	  IsSet[i][j]=1;
	  ij2fzfz(i, j, &newa1, &newa2, &newb1, &newb2, N1, N2);
	  if(newa1!=a1 || newa2!=a2 || newb1!=b1 || newb2!=b2) error("test_fzij: function ij2fzfz not the inverse of fzfz2ij."); 
	}
  for(i=0; i<N1*N2; i++) for(j=0; j<N1*N2; j++) if (IsSet[i][j]!=1) error("Matrix [i,j] not fully populated.");

  free_matrix(IsSet,0,N1*N2-1,0,N1*N2-1);
}
