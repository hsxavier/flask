#include "Cosmology.hpp"
#include "Utilities.hpp"
#include "Integral.hpp"
#include "interpol.h"
#include "ParameterList.hpp"
#include <math.h>

/*** Methods of the Cosmology class ***/

Cosmology::Cosmology() {
  curv=1;
  //curv=0;       // Cosmological distances in ComDist function are calculated assuming flat universe.
  lowz=1;       // Cosmological distances in ComDist function are calculated ignoring radiation contributions.
  c=299792.458; /* km/s     */
  H100=100;     /* km/s/Mpc */
  NOTSETd=-666;
  NOTSETi=-666; 

  Om      = NOTSETd;    
  Ode     = NOTSETd;
  Ob      = NOTSETd;
  Onu     = NOTSETd;
  Omh2    = NOTSETd;    
  Odeh2   = NOTSETd;
  Obh2    = NOTSETd;
  Onuh2   = NOTSETd;
  Ok      = NOTSETd;
  wde     = NOTSETd;
  h       = NOTSETd;
  Nnu     = NOTSETi;
  deltaH2 = NOTSETd;
  ns      = NOTSETd;
  galdens = NOTSETd;
}

void Cosmology::load(ParameterList *config) {
  using namespace ParDef;
  //Omh2    = config->readd("OMEGA_mh2");
  Om      = config->readd("OMEGA_m");
  //Odeh2   = config->readd("OMEGA_Lh2");
  Ode     = config->readd("OMEGA_L");
  wde     = config->readd("W_de");
  //Obh2    = config->readd("OMEGA_bh2");
  //Onuh2   = config->readd("OMEGA_nuh2");
  //Nnu     = config->readi("N_nu");
  //h       = config->readd("HUBBLE");
  //deltaH2 = config->readd("PK_DELTAH")*config->readd("PK_DELTAH");
  //ns      = config->readd("PK_NS");
  galdens   = config->readd("GALDENSITY");
  SetOther();
}

void Cosmology::SetOther() {
  //if (h==NOTSETd)       error("Cosmology::SetOther: Hubble constant must be set.");
  if (wde==NOTSETd)     error("Cosmology::SetOther: equation of state wde must be set.");
  //if (Nnu==NOTSETi)     error("Cosmology::SetOther: number of massive neutrinos Nnu must be set.");
  //if (deltaH2==NOTSETi) error("Cosmology::SetOther: power spectrum normalization deltaH2 must be set.");
  //if (ns==NOTSETi)      error("Cosmology::SetOther: power spectrum index ns must be set.");
  /*
  // Set Omega_matter:
  if (Om==NOTSETd && Omh2==NOTSETd) error("Cosmology::SetOther: Om or Omh2 must be set.");
  if (Om!=NOTSETd && Omh2!=NOTSETd && Om!=Omh2/h/h) error("Cosmology::SetOther: conflicting values of Om and Omh2.");
  if (Om==NOTSETd) Om   = Omh2/h/h;
  else             Omh2 = Om  *h*h;
  // Set Omega_DE:
  if (Ode==NOTSETd && Odeh2==NOTSETd) error("Cosmology::SetOther: Ode or Odeh2 must be set.");
  if (Ode!=NOTSETd && Odeh2!=NOTSETd && Ode!=Odeh2/h/h) error("Cosmology::SetOther: conflicting values of Ode and Odeh2.");
  if (Ode==NOTSETd) Ode   = Odeh2/h/h;
  else              Odeh2 = Ode  *h*h;
  // Set Omega_barion:
  if (Ob==NOTSETd && Obh2==NOTSETd) error("Cosmology::SetOther: Ob or Obh2 must be set.");
  if (Ob!=NOTSETd && Obh2!=NOTSETd && Ob!=Obh2/h/h) error("Cosmology::SetOther: conflicting values of Ob and Obh2.");
  if (Ob==NOTSETd) Ob   = Obh2/h/h;
  else             Obh2 = Ob  *h*h;
  // Set Omega_neutrino:
  if (Onu==NOTSETd && Onuh2==NOTSETd) error("Cosmology::SetOther: Onu or Onuh2 must be set.");
  if (Onu!=NOTSETd && Onuh2!=NOTSETd && Onu!=Onuh2/h/h) error("Cosmology::SetOther: conflicting values of Onu and Onuh2.");
  if (Onu==NOTSETd) Onu   = Onuh2/h/h;
  else              Onuh2 = Onu  *h*h;
  */
  // Set Omega_curvature:
  if (Om!=NOTSETd && Ode!=NOTSETd) Ok = 1.0-Om-Ode;
  else error("Cosmology::SetOther: cannot set curvature due to missing Om and/or Ode.");
  // Check spatial geometry:
  if (fabs(Om+Ode-1.0)>0.0001 && curv==0) {
    char message[100];
    sprintf(message, "Cosmology::SetOther: ComDist assumes flat universe but Om+Ode= %g.", Om+Ode);
    warning(message);
  }
}

void Cosmology::show (std::ostream * output) {
  if (output!=&std::cout) warning("Cosmology::show: output to file not implemented.");

  printf("%10s %10g %50s\n", "Om:",Om,"DM+barions+massive neutrinos density parameter.");
  printf("%10s %10g %50s\n", "Ob:",Ob,"Barion density parameter.");
  printf("%10s %10.3g %50s\n", "Onu:",Onu,"Massive neutrino density parameter.");
  printf("%10s %10g %50s\n", "Ode:",Ode,"Dark energy density parameter.");
  printf("%10s %10g %50s\n", "Omh2:",Omh2,"Total matter density parameter times h^2.");
  printf("%10s %10g %50s\n", "Obh2:",Obh2,"Barion density parameter times h^2.");
  printf("%10s %10g %50s\n", "Onuh2:",Onuh2,"Massive neutrino density parameter times h^2.");
  printf("%10s %10g %50s\n", "Odeh2:",Odeh2,"Dark energy density parameter times h^2.");
  printf("%10s %10g %50s\n", "Ok:",Ok,"Curvature density parameter.");
  printf("%10s %10g %50s\n", "wde:",wde,"Dark energy equation of state.");
  printf("%10s %10d %50s\n", "Nnu:",Nnu,"Number of massive species of neutrinos.");
  printf("%10s %10g %50s\n", "deltaH2:",deltaH2,"Power spectrum normalization.");
  printf("%10s %10g %50s\n", "ns:",ns,"Power spectrum index.");
  printf("%10s %10g %50s\n", "Tcmb:", 2.728, "Hard-wired in WayneHuPowerClass.cpp.");
  printf("%10s %10g %50s\n", "h:",h,"Hubble constant in units of 100km/s/Mpc.");
  printf("%10s %10g %50s\n", "H100:",H100,"100 km/s/Mpc.");
  printf("%10s %10g %50s\n", "c:",c,"Speed of light in km/s.");
  printf("%10s %10d %50s\n", "curv:",curv,"Hard-wired, safe-guard DistCom function.");
  printf("%10s %10d %50s\n", "lowz:",lowz,"Hard-wired, no implications.");
  printf("%10s %10g %50s\n", "galdens:",galdens,"Galaxy comoving number density in (h^-1 Mpc)^-3.");
}


/*** Separate cosmological functions ***/

double Eh(Cosmology *p, double z) {
  return sqrt(p->Om*pow(1+z,3) + p->Ok*pow(1+z,2) + p->Ode*pow(1+z,3*(1+p->wde)));
}
double ComDistIntegrand(Cosmology *p, double z) {
  return 1.0/Eh(p, z);
}
// Returns the comoving distance along the line of sight in h^-1 Mpc:
double ComDist(Cosmology *p, double z) {
  const int ngrid=500;
  const double zmin=0.0, zmax=8.0, maxerr=0.1;
  static bool init=0;
  static double zgrid[ngrid+1], dgrid[ngrid+1];
  int i;
  double dist;
  // Initialize:
  if (init==0) {
    for (i=0; i<=ngrid; i++) {
      zgrid[i] = zmin + i*((zmax-zmin)/ngrid);
      dgrid[i] = p->c/p->H100*qromb(ComDistIntegrand, 0.0, zgrid[i], p);
    }
    init=1;
  }
  // Return from interpolation:
  dist=Interpol(zgrid, ngrid+1, dgrid, z);
  return dist;  
}

// Derivative of Comoving distance with respect to redshift:
double dXidz(Cosmology *p, double z) {
  return p->c/p->H100*ComDistIntegrand(p, z);
}

// Redshift selection function for projection:
double zSelection(double z, double z0) {
  return 1.0;
}
double ProjDensityIntegrand(double z, double z0, Cosmology *p) {
  return zSelection(z,z0) * pow(ComDist(p,z),2) * dXidz(p,z);
}
double ProjDensityIntegrand(double z, Cosmology *p) {
  return pow(ComDist(p,z),2) * dXidz(p,z);
}
double ProjDensity(double z0, double zmin, double zmax, Cosmology *p) {
  return p->galdens * qromb(ProjDensityIntegrand, zmin, zmax, z0, p);
}

// It is the modulus squared of the fourier transform of a 3D Tophat function with normalization 1. 
double TophatWk2(double kR) {
  if (kR<0.001) return 1.0; // Error caused by this approximation should be less than 1e-6.
  return 9.0*pow(sin(kR)-kR*cos(kR), 2)/pow(kR,6);
}
