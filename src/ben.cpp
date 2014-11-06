#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Utilities.hpp"
#include "Utilities.cpp"

#define PI M_PI
#define LN2 0.69314718

int main() {
  double **cl;
  long ncl, status;

  cl = LoadTable<double>("cltest.dat", &ncl, &status);
  
  return 0;
}

void bisect_interpol(double *a,int dim,double x,int *erg)
{
  int l=0,u=dim-1,m;

  if ((x<a[l]) || (x>a[u])) {
    printf("Error: index out of range: %g not in %g-%g\n",x,a[l],a[u]);
    exit(-1);
  }
  while (u-l>1) {
    m=(u+l)/2;
    if (x>a[m]) l=m;
    else u=m;
  }
  *erg=l;
}

double interpol_linear_extra_1D(double *ax,int NX,double *z,double x)
{
  if (x<=ax[0]) return((z[1]-z[0])/(ax[1]-ax[0])*(x-ax[0])+z[0]);
  else if (x>=ax[NX-1]) return((z[NX-1]-z[NX-2])/(ax[NX-1]-ax[NX-2])*(x-ax[NX-1])+z[NX-1]);
  else {
    int ix;
    bisect_interpol(ax,NX,x,&ix);
    double t=(x-ax[ix])/(ax[ix+1]-ax[ix]);
    return((1-t)*z[ix]+t*z[ix+1]);
  }
}

double *bj_alloc(int dim) // Allocate an array [0...dim-1] of doubles (zero-initialized).
{
  double *arr=(double*)calloc(dim,sizeof(double));
  return arr;
}

void math_cdgamma(fftw_complex x, fftw_complex *res)
{
 double xr, xi, wr, wi, ur, ui, vr, vi, yr, yi, t;

 xr = (double) x[0];
 xi = (double) x[1];
 if (xr < 0) {
   wr = 1 - xr;
   wi = -xi;
 } 
 else {
   wr = xr;
   wi = xi;
 }
 ur = wr + 6.00009857740312429;
 vr = ur * (wr + 4.99999857982434025) - wi * wi;
 vi = wi * (wr + 4.99999857982434025) + ur * wi;
 yr = ur * 13.2280130755055088 + vr * 66.2756400966213521 + 0.293729529320536228;
 yi = wi * 13.2280130755055088 + vi * 66.2756400966213521;
 ur = vr * (wr + 4.00000003016801681) - vi * wi;
 ui = vi * (wr + 4.00000003016801681) + vr * wi;
 vr = ur * (wr + 2.99999999944915534) - ui * wi;
 vi = ui * (wr + 2.99999999944915534) + ur * wi;
 yr += ur * 91.1395751189899762 + vr * 47.3821439163096063;
 yi += ui * 91.1395751189899762 + vi * 47.3821439163096063;
 ur = vr * (wr + 2.00000000000603851) - vi * wi;
 ui = vi * (wr + 2.00000000000603851) + vr * wi;
 vr = ur * (wr + 0.999999999999975753) - ui * wi;
 vi = ui * (wr + 0.999999999999975753) + ur * wi;
 yr += ur * 10.5400280458730808 + vr;
 yi += ui * 10.5400280458730808 + vi;
 ur = vr * wr - vi * wi;
 ui = vi * wr + vr * wi;
 t = ur * ur + ui * ui;
 vr = yr * ur + yi * ui + t * 0.0327673720261526849;
 vi = yi * ur - yr * ui;
 yr = wr + 7.31790632447016203;
 ur = log(yr * yr + wi * wi) * 0.5 - 1;
 ui = atan2(wi, yr);
 yr = exp(ur * (wr - 0.5) - ui * wi - 3.48064577727581257) / t;
 yi = ui * (wr - 0.5) + ur * wi;
 ur = yr * cos(yi);
 ui = yr * sin(yi);
 yr = ur * vr - ui * vi;
 yi = ui * vr + ur * vi;
 if (xr < 0) {
   wr = xr * 3.14159265358979324;
   wi = exp(xi * 3.14159265358979324);
   vi = 1 / wi;
   ur = (vi + wi) * sin(wr);
   ui = (vi - wi) * cos(wr);
   vr = ur * yr + ui * yi;
   vi = ui * yr - ur * yi;
   ur = 6.2831853071795862 / (vr * vr + vi * vi);
   yr = ur * vr;
   yi = ur * vi;
 }
 (*res)[0]=yr; 
 (*res)[1]=yi;
}

// HX:                angle?     kernel (complex num)     [0,0]      2
//Convolution kernel for Hankel-Transform with bias q, Bessel_mu
void Hankel_Kernel_FT(double x, fftw_complex *res, double *arg, int argc)
{
  //arg[0]: q, arg[1]:mu
  fftw_complex a1,a2,g1,g2;
  int mu;
  double mod,xln2,si,co,d1,d2,pref,q;

  q=arg[0];
  mu=(int) (arg[1]+0.1);

  //arguments for complex cosmopar.Gamma
  a1[0]=0.5*(1.0+mu+q);
  a2[0]=0.5*(1.0+mu-q);
  a1[1]=0.5*x; a2[1]=-a1[1];

  math_cdgamma(a1,&g1); 
  math_cdgamma(a2,&g2);

  xln2=x*LN2;
  si=sin(xln2); co=cos(xln2);
  d1=g1[0]*g2[0]+g1[1]*g2[1]; //Re
  d2=g1[1]*g2[0]-g1[0]*g2[1]; //Im
  mod=g2[0]*g2[0]+g2[1]*g2[1];
  pref=exp(LN2*q)/mod;
  
  (*res)[0]=pref*(co*d1-si*d2);
  (*res)[1]=pref*(si*d1+co*d2);
}

// HX: Hankel transf. f(x) out,     x out,   N_x out,       f(x) in,          x in,     N_x in,      x in min,        x in max.  
void hankeltrafo(double *xi,double *theta,int N_OUT,double *logsignal,double *logarg,int N_IN,double samplemin,double samplemax)
{
  double dlnl,loglmax,loglmin,lnell,kk,arg[2],lnrc;
  double *lP;
  fftw_plan plan1,plan;
  fftw_complex *f_lP,*conv;
  fftw_complex kernel;
  int i,nc;

  lP=(double*)fftw_malloc(N_OUT*sizeof(double));               // HX: ell*P sampled in log(ell) bins
  f_lP=(fftw_complex*)fftw_malloc((N_OUT/2+1)*sizeof(fftw_complex)); // HX: Fourier transform of the above.
  conv=(fftw_complex*)fftw_malloc((N_OUT/2+1)*sizeof(fftw_complex)); 

  // HX:                          in  out                          
  plan=fftw_plan_dft_r2c_1d(N_OUT,lP,f_lP,FFTW_ESTIMATE);
  // HX:                           in  out
  plan1=fftw_plan_dft_c2r_1d(N_OUT,conv,lP,FFTW_ESTIMATE);

  loglmax=samplemax;
  loglmin=samplemin;   // add margin due to strange nc convention
  dlnl=(loglmax-loglmin)/(1.0*N_OUT-1);
  lnrc=0.5*(loglmax+loglmin);
  nc=N_OUT/2+1;

  //Power spectrum on logarithmic bins
  for(i=0;i<N_OUT;i++) {
    lnell=lnrc+(i-nc)*dlnl;
    if (lnell<=logarg[0]) lP[i]=exp(lnell+logsignal[0]+(logsignal[1]-logsignal[0])/(logarg[1]-logarg[0])*(lnell-logarg[0]));
    else if (lnell>=logarg[N_IN-1]) lP[i]=exp(lnell+logsignal[N_IN-2]+(logsignal[N_IN-1]-logsignal[N_IN-2])/(logarg[N_IN-1]-logarg[N_IN-2])*(lnell-logarg[N_IN-2]));
    else lP[i]=exp(lnell+interpol_linear_extra_1D(logarg,N_IN,logsignal,lnell));   //sample ell*P on logarithmic bins
  }

  //Go to log-Fourier-space
  fftw_execute(plan); // HX: It seems to me that I am going to correlation function.

  arg[0]=0; // unused
  arg[1]=0; // J_0

  // HX: What is this?
  //perform the convolution, negative sign for kernel (cc!)
  for(i=0;i<N_OUT/2+1;i++) {
    kk=2*PI*i/(dlnl*N_OUT);              // HX: seems like an angle value?
    Hankel_Kernel_FT(kk,&kernel,arg,2);  // HX: Compute the factor kernel(kk).
    conv[i][0]=f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1]; // HX: Re part of f_lP*kernel. Parei aqui.
    conv[i][1]=f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1]; // HX: Im part of f_lP*kernel.
  }
  conv[0][1]=0; //Nyquist- and 0-frequency-components are double!
  conv[N_OUT/2][1]=0;

  //go back to double space, i labels log-bins in theta
  fftw_execute(plan1); // HX: Go back to the Power-spectrum?

  for(i=0;i<N_OUT;i++) {
    theta[N_OUT-i-1]=exp((nc-i)*dlnl-lnrc);
    xi[N_OUT-i-1]=lP[i]/(theta[N_OUT-i-1]*2*PI*N_OUT) ;   //sample ell*P_2 on logarithmic bins
  }

  //Clean up
  fftw_free(conv);
  fftw_free(lP);
  fftw_free(f_lP);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
  return;
}



// takes & returns logarithmic ell and ps (in-place trafo)
void convertpstolognormal(double *ell,double *ps,int NBIN,double shift)
{
  int i,flag=0;
  const int NT=NBIN;
  int NT2=NBIN;

  double limlow,limhig,vallow,valhig;
  double *theta=bj_alloc(NT);
  double *xi=bj_alloc(NT);

  vallow=ps[0]+ell[0];            // find FFT limits (ps & ell logged)
  valhig=ps[NBIN-1]+ell[NBIN-1];  // HX: actually ln(ps)+ln(ell)=ln(ell*ps)
  i=NBIN-1;
  while (valhig-vallow<-7.0) {    // HX: -7.0 should be #defined or something.
    i--;
    valhig=ps[i]+ell[i];
    if (i==0) {
      printf("Error in 'convertpstolognormal': Could not find FFT limits!\n");
      exit(-1);
    }
  }
  limhig=ell[i];
  i=0;
  while (valhig-vallow>7.0) {     // HX: 7.0 should be #defined of something.
    i++;                          // HX: It is not clear to me why this is done and what is done. 
    vallow=ps[i]+ell[i]; 
    if (i==NBIN-1) {
      printf("Error in 'convertpstolognormal': Could not find FFT limits!\n");
      exit(-1);
    }
  }
  limlow=ell[i];
  printf("Limits 1. FFT: %g  %g\n",exp(limlow),exp(limhig));

  hankeltrafo(xi,theta,NT,ps,ell,NBIN,limlow,limhig); // HX: This seems return the correlation function.

  for(i=0;i<NT;i++) {
    if (xi[i]<=0.0) {
      xi[i]=1.e-10;  // avoid negative xi (can be avoided by modifying hankeltrafo routine)
      if (!flag) printf("Encountered negative xi at theta=%8.2f deg\n",theta[i]/PI*180.);
      flag=1;
    }
    /*
    if (xi[i]<0.0) {
      NT2=i;
      printf("Encountered negative xi at theta=%.2g deg -> limit data to %i bins!\n",theta[i]/PI*180.,NT2);
      break;
    }
    */
    //printf("%g  %g\n",theta[i],xi[i]); // HX: Get the gaussian correlation function.
    xi[i]=log(1.+xi[i]/(shift*shift));  
    xi[i]=log(xi[i]);         //transform to log-space
    theta[i]=log(theta[i]);   //transform to log-space
  }
  limlow=theta[0];
  limhig=theta[NT2-1];
  printf("Limits 2. FFT: %g  %g\n",exp(limlow),exp(limhig));

  hankeltrafo(ps,ell,NBIN,xi,theta,NT2,limlow,limhig); // Get the gaussian Power-spectrum.

  for(i=0;i<NBIN;i++) {
    ps[i]*=(4.*PI*PI);   // correct for asymmetric Fourier trafo
    ps[i]=log(ps[i]);    //transform to log-space
    ell[i]=log(ell[i]);  //transform to log-space
  }

  free(theta);
  free(xi);
  return;
}
