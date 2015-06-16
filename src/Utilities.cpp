/***************************************************************
2012-09-18: Utility functions: error handling. 
***************************************************************/

#include "Utilities.hpp"
#include <iostream>
#include <cstdlib>   // to use exit()
#include <iomanip>   // for setw()
#include <ctime>

/****** Replace substring 'from' to 'to' in string 'str' ******/
bool StrReplace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos) return false;
    str.replace(start_pos, from.length(), to);
    return true;
}


/****** Print message and measure time elapsed from message printing to next function call with empty parameters ******/
void Announce(std::string message) {
  static time_t start;
  int space=58;
  double diff;

  // Announce the beginning of a process:
  if (message!="done") {
    //std::cout.width(space);
    std::cout << std::left << std::setw(space) << message << std::right; std::cout.flush();
    start = time(NULL);
  }

  // Announce the end of a process: 
  else {
    diff = difftime(time(NULL), start);
    printf("done.  (%gs)\n", diff);
  }
}

/****** Import vertical vectors from file ******/
// No memory allocation is done here. It must be done before.  
void ImportVecs(double **matriz, long length, long nvecs, const char *filename) {
   
   /* Declaração das variáveis */
   FILE *arq;
   char message[100];
   long i, j;
   
   if ((arq=fopen(filename, "r"))==NULL) {
     sprintf(message,"ImportVecs: cannot open '%s' file.", filename);
     error(message);
   }
    
   /* Leitura da matriz */
   for (i=0; i<length; i++)
     for (j=0; j<nvecs; j++)
       fscanf(arq, "%lf", &matriz[j][i]);
   
   fclose(arq);
}


// Error handling functions:
void error (const std::string message) {
  //using std::cout;
  using std::cerr;
  //cout << "\n!! ERROR!  " << message << " !!\n";
  cerr << "\n!! ERROR!  " << message << " !!\n";
  exit(1);
}
int warning (const std::string message) {
  using std::cout;
  using std::cerr;
  static int counter=0;
  if (message!="count") {
    //cout << "\n!! WARNING!  " << message << " !!\n";
    cerr << "\n!! WARNING!  " << message << " !!\n";
    counter++;
    return 0;
  }
  return counter;
}


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.0e-12        // Altered to handle double precision
#define RNMX (1.0-EPS)
// Random generator function from numerical recipes originally called 'ran2':
double random(long *idum) {
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#include <math.h>
// Gaussian random number generator from Numerical Recipes
double gasdev(long *idum) {
  double random(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*random(idum)-1.0;
      v2=2.0*random(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software #?w,(1. */

/*** Pads a string containing a number with zeroes on the left ***/
std::string ZeroPad(int num, int max) {
  std::stringstream ss;
  int ndigits=1;
  
  while (max >= 10) {max=max/10; ndigits++;}
  
  
  ss << std::setfill('0') << std::setw(ndigits) << num;
  return ss.str();
}
