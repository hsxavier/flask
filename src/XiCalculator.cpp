#include <iostream>
#include "Utilities.hpp"
#include <cmath>
#include <omp.h>                // For OpenMP functions, not pragmas.

int main (int argc, char *argv[]) {
  const double Degrees=M_PI/180.0;
  const int RAcol=1, DECcol=0;
  const double MinAng=0, MaxAng=M_PI;
  double DeltaAng;
  int MaxThreads, p;
  std::string infile, outfile;
  std::ofstream output;
  long Nrows, Ncols, i, j, bin, nXi, k, kmax;
  double **cat, **Xi, *fXi, *Ang, dist, *XiCounts, *wrapper[3];
  int **Counts, VARcol, *fCounts;

  // Get input from command line:
  if (argc<=4) {
    printf("USAGE: XiCalculator <INPUT FILE> <FIELD COL> <DELTA_ANG (DEG)> <OUTPUT FILE>\n");
    printf("                                  1          2               N              \n");
    printf("<INPUT FILE> columns must be: DEC [deg], RA [deg], ..., FIELD VALUE, ...    \n");
    return 0;
  }
  infile.assign(argv[1]);
  sscanf(argv[2],"%d",  &VARcol);
  sscanf(argv[3],"%lf", &DeltaAng);
  outfile.assign(argv[4]);
  printf("                                            [0] [1]        [%d]\n", VARcol);
  printf("Will compute Xi with resolution %g deg for (DEC, RA, ..., FIELD) in file %s:\n", 
	 DeltaAng, infile.c_str());
  
  // Initialize bunch of things:
  VARcol   = VARcol-1;
  DeltaAng = DeltaAng*Degrees;
  std::cout << "Loading catalog...\n";
  cat      = LoadTable<double>(infile, &Nrows, &Ncols, 0, 1);
  std::cout << "Preparing for Xi calculations...\n";
  MaxThreads = omp_get_max_threads();
  nXi        = (long)((MaxAng-MinAng)/DeltaAng);
  Xi         = matrix<double>(0, MaxThreads-1, 0, nXi);
  Counts     = matrix<int>   (0, MaxThreads-1, 0, nXi);
#pragma omp parallel for private(i)
  for(p=0; p<MaxThreads; p++) {
    for (i=0; i<=nXi; i++) {
      Xi[p][i]     = 0.0;
      Counts[p][i] = 0;
    }
  }
  
  // Go from degrees to radians:
#pragma omp parallel for
  for(i=0; i<Nrows; i++) {
    cat[i][ RAcol] = cat[i][ RAcol]*Degrees;
    cat[i][DECcol] = cat[i][DECcol]*Degrees;
  }



  // LOOP over pairs, summing over cross multiplications:
  std::cout << "Calculating Xi...\n";
  kmax = (Nrows*(Nrows-1))/2;
#pragma omp parallel for schedule(static) private(p, i, j, dist, bin)
  for (k=0; k<kmax; k++) {
    p = omp_get_thread_num();
    i = (int)((sqrt(8.0*k+1.0)-1.0)/2.0) + 1;
    j = k-((i-1)*i)/2;
     
    dist = acos( sin(cat[i][DECcol])*sin(cat[j][DECcol]) + cos(cat[i][DECcol])*cos(cat[j][DECcol])*cos(cat[i][RAcol]-cat[j][RAcol]) );
    bin  = (long)((dist - MinAng)/DeltaAng);
    Xi[p][bin] += cat[i][VARcol]*cat[j][VARcol];
    Counts[p][bin]++;
  }
  free_matrix(cat, 0, Nrows-1, 0, Ncols-1);

  std::cout << "Post processing...\n";
  // Initialize new vectors:
  fXi      = vector<double>(0, nXi);
  fCounts  = vector<int>   (0, nXi);
  Ang      = vector<double>(0, nXi);
  XiCounts = vector<double>(0, nXi);
  for (i=0; i<=nXi; i++) {
    fXi[i]     = 0.0;
    fCounts[i] = 0;
  }
  // Take total over partial sum:
#pragma omp parallel for private(p)
  for (i=0; i<=nXi; i++) {
    for(p=0; p<MaxThreads; p++) {  
      fXi[i]     += Xi[p][i];
      fCounts[i] += Counts[p][i];
    }
  }
  // Take the mean:
  for (i=0; i<=nXi; i++) {
    if (fCounts[i]>0) fXi[i] = fXi[i]/((double)fCounts[i]); 
    Ang[i]      = MinAng + ((double)i+0.5)*DeltaAng;
    XiCounts[i] = fCounts[i] / sin(Ang[i]); 
  }
  if (fXi[nXi]>0) warning("fXi[nXi]>0 is weird.");

  // Export Xi:
  wrapper[0]=Ang;
  wrapper[1]=fXi;
  output.open(outfile.c_str());
  if (!output.is_open()) error("Cannot use output file "+outfile);
  PrintVecs(wrapper, nXi, 2, &output);
  output.close();

  // Free memory and exit:
  free_matrix(Xi,       0, MaxThreads-1, 0, nXi);
  free_matrix(Counts,   0, MaxThreads-1, 0, nXi);
  free_vector(XiCounts, 0, nXi);
  free_vector(Ang,      0, nXi);
  free_vector(fXi,      0, nXi);
  free_vector(fCounts,  0, nXi);
  return 0;
}
