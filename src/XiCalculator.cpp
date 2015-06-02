#include <iostream>
#include "Utilities.hpp"
#include <cmath>
#include <omp.h> // For OpenMP functions, not pragmas.


int main (int argc, char *argv[]) {
  const double Degrees=M_PI/180.0;
  const int RAcol=0, DECcol=1;
  const double MinAng=0, MaxAng=M_PI;
  std::string datafile, randfile, outfile;
  std::ofstream output;
  double DeltaAng;
  int MaxThreads, p, VARcol;
  long dataNrows, dataNcols, randNrows, randNcols, i, j, bin, nXi, k, kmax;
  double **data, **rand, **PartFF, *FF, *Ang, dist, *wrapper[3], *LandySzalay, factor1, factor2;
  int **PartDD, *DD, **PartDR, *DR, **PartRR, *RR;


  /****************************************/
  /*** Reading input and preparing data ***/
  /****************************************/

  // Get input from command line:
  if (argc<=5) {
    printf("USAGE: XiCalculator <DATA FILE> <FIELD COL> <RAND FILE> <DELTA_ANG (DEG)> <OUTPUT FILE>\n");
    printf("                                 %d        %d               N                          \n", RAcol, DECcol);
    printf("<DATA FILE> columns must be: RA [deg], DEC [deg], ..., FIELD VALUE, ...                \n");
    printf("<RAND FILE> columns must be: RA [deg], DEC [deg].                                      \n");
    return 0;
  }
  datafile.assign(argv[1]);
  sscanf(argv[2],"%d",  &VARcol);
  randfile.assign(argv[3]);
  sscanf(argv[4],"%lf", &DeltaAng);
  outfile.assign(argv[5]);
  printf("                                             [%d] [%d]       [%d]\n", DECcol, RAcol, VARcol);
  printf("Will compute Xi with resolution %g deg for (DEC, RA, ..., FIELD) in file %s using %s as random sample:\n", 
	 DeltaAng, datafile.c_str(), randfile.c_str());
  
  // Initialize bunch of things and load files:
  std::cout << "Preparing for calculations...\n";
  VARcol     = VARcol-1;
  DeltaAng   = DeltaAng*Degrees;
  MaxThreads = omp_get_max_threads();
  nXi        = (long)((MaxAng-MinAng)/DeltaAng);
  data       = LoadTable<double>(datafile, &dataNrows, &dataNcols, 0, 1);
  rand       = LoadTable<double>(randfile, &randNrows, &randNcols, 0, 1);
  
  // Go from degrees to radians:
#pragma omp parallel for
  for(i=0; i<dataNrows; i++) {
    data[i][ RAcol] = data[i][ RAcol]*Degrees;
    data[i][DECcol] = data[i][DECcol]*Degrees;
  }
#pragma omp parallel for
  for(i=0; i<randNrows; i++) {
    rand[i][ RAcol] = rand[i][ RAcol]*Degrees;
    rand[i][DECcol] = rand[i][DECcol]*Degrees;
  }


  /******************************************************/
  /*** Compute Data pair counts and Field correlation ***/
  /******************************************************/

  // Initialize arrays for partial sums:
  PartFF     = matrix<double>(0, MaxThreads-1, 0, nXi);
  PartDD     = matrix<int>   (0, MaxThreads-1, 0, nXi);
#pragma omp parallel for private(i)
  for(p=0; p<MaxThreads; p++) {
    for (i=0; i<=nXi; i++) {
      PartFF[p][i] = 0.0;
      PartDD[p][i] = 0;
    }
  }

  // LOOP over DATA pairs (DD and FF), summing over cross multiplications of F:
  std::cout << "Summing over pairs DD and FF...\n";
  kmax = (dataNrows*(dataNrows-1))/2;
#pragma omp parallel for schedule(static) private(p, i, j, dist, bin)
  for (k=0; k<kmax; k++) {
    p = omp_get_thread_num();
    i = (int)((sqrt(8.0*k+1.0)-1.0)/2.0) + 1;
    j = k-((i-1)*i)/2;
    
    // Find angular distance bin:
    dist = acos( sin(data[i][DECcol])*sin(data[j][DECcol]) + cos(data[i][DECcol])*cos(data[j][DECcol])*cos(data[i][RAcol]-data[j][RAcol]) );
    bin  = (long)((dist - MinAng)/DeltaAng);
    // Add value to the bin:
    PartFF[p][bin] += data[i][VARcol]*data[j][VARcol];
    PartDD[p][bin] ++;
  }

  // Prepare to sum over partial sums:
  std::cout << "Finalizing summation...\n";
  FF       = vector<double>(0, nXi);
  DD       = vector<int>   (0, nXi);
  for (i=0; i<=nXi; i++) {
    FF[i] = 0.0;
    DD[i] = 0;
  }
  // Total over partial sums:
#pragma omp parallel for private(p)
  for (i=0; i<=nXi; i++) {
    for(p=0; p<MaxThreads; p++) {  
      FF[i] += PartFF[p][i];
      DD[i] += PartDD[p][i];
    }
  }
  free_matrix(PartFF, 0, MaxThreads-1, 0, nXi);
  free_matrix(PartDD, 0, MaxThreads-1, 0, nXi);  
  // Normalize the sum FF by the number of points:
  for (i=0; i<=nXi; i++) if (DD[i]>0) FF[i] = FF[i]/((double)DD[i]); 
  if (FF[nXi]>0) warning("FF[nXi]>0 is weird.");


  /**********************************/
  /*** Compute Random pair counts ***/
  /**********************************/

  // Initialize arrays for partial sums:
  PartRR = matrix<int>(0, MaxThreads-1, 0, nXi);
#pragma omp parallel for private(i)
  for(p=0; p<MaxThreads; p++) {
    for (i=0; i<=nXi; i++) PartRR[p][i] = 0;
  }
  
  // LOOP over RAND pairs (RR):
  std::cout << "Summing over pairs RR...\n";
  kmax = (randNrows*(randNrows-1))/2;
#pragma omp parallel for schedule(static) private(p, i, j, dist, bin)
  for (k=0; k<kmax; k++) {
    p = omp_get_thread_num();
    i = (int)((sqrt(8.0*k+1.0)-1.0)/2.0) + 1;
    j = k-((i-1)*i)/2;
    
    // Find angular distance bin:
    dist = acos(sin(rand[i][DECcol])*sin(rand[j][DECcol]) + cos(rand[i][DECcol])*cos(rand[j][DECcol])*cos(rand[i][RAcol]-rand[j][RAcol]));
    bin  = (long)((dist - MinAng)/DeltaAng);
    // Add value to the bin:
    PartRR[p][bin]++;
  }

  // Prepare to sum over partial sums:
  std::cout << "Finalizing summation...\n";
  RR = vector<int>(0, nXi);
  for (i=0; i<=nXi; i++) RR[i] = 0;
  // Total over partial sums:
#pragma omp parallel for private(p)
  for (i=0; i<=nXi; i++) {
    for(p=0; p<MaxThreads; p++) RR[i] += PartRR[p][i];
  }
  free_matrix(PartRR, 0, MaxThreads-1, 0, nXi);  


  /***************************************/
  /*** Compute Data-Random pair counts ***/
  /***************************************/
  
  // Initialize arrays for partial sums:
  PartDR = matrix<int>(0, MaxThreads-1, 0, nXi);
#pragma omp parallel for private(i)
  for(p=0; p<MaxThreads; p++) {
    for (i=0; i<=nXi; i++) PartDR[p][i] = 0;
  }
  
  // LOOP over data-random pairs (DR):
  std::cout << "Summing over pairs DR...\n";
  kmax = dataNrows*randNrows;
#pragma omp parallel for schedule(static) private(p, i, j, dist, bin)
  for (k=0; k<kmax; k++) {
    p = omp_get_thread_num();
    i = k/randNrows; // D index
    j = k%randNrows; // R index
    
    // Find angular distance bin:
    dist = acos(sin(data[i][DECcol])*sin(rand[j][DECcol]) + cos(data[i][DECcol])*cos(rand[j][DECcol])*cos(data[i][RAcol]-rand[j][RAcol]));
    bin  = (long)((dist - MinAng)/DeltaAng);
    // Add value to the bin:
    PartDR[p][bin]++;
  }

  // Prepare to sum over partial sums:
  std::cout << "Finalizing summation...\n";
  DR = vector<int>(0, nXi);
  for (i=0; i<=nXi; i++) DR[i] = 0;
  // Total over partial sums:
#pragma omp parallel for private(p)
  for (i=0; i<=nXi; i++) {
    for(p=0; p<MaxThreads; p++) DR[i] += PartDR[p][i];
  }
  free_matrix(PartDR, 0, MaxThreads-1, 0, nXi);  


  /**************************************************************/
  /*** Compute density correlation function with Landy-Szalay ***/
  /**************************************************************/
  
  printf("Computing Landy-Szalay correlation estimator...\n");
  LandySzalay = vector<double>(0, nXi);
  factor1     = ((double)(randNrows*(randNrows-1))) / ((double)(dataNrows*(dataNrows-1)));
  factor2     = ((double)           (randNrows-1))  / ((double) dataNrows               ); 
  for (i=0; i<nXi; i++) {
    if (RR[i]>0) LandySzalay[i] = factor1*((double)DD[i])/((double)RR[i]) - factor2*((double)DR[i])/((double)RR[i]) + 1.0; 
  }
  
  
  /********************************************/
  /*** Clear data and rand catalog memories ***/
  /********************************************/

  free_matrix(data, 0, dataNrows-1, 0, dataNcols-1);
  free_matrix(rand, 0, randNrows-1, 0, randNcols-1);  


  /*************************/
  /*** Exporting results ***/
  /*************************/

  // Generate list of separating angles:
  printf("Outputting results...\n");
  Ang = vector<double>(0, nXi);
  for (i=0; i<=nXi; i++) Ang[i] = (MinAng + ((double)i+0.5)*DeltaAng)/Degrees;

  // Export correlation functions:
  wrapper[0] = Ang;
  wrapper[1] = LandySzalay;
  wrapper[2] = FF;
  output.open(outfile.c_str());
  if (!output.is_open()) error("Cannot use output file "+outfile);
  PrintVecs(wrapper, nXi, 3, &output);
  output.close();

  // Free memory and exit:
  free_vector(Ang,    0, nXi);
  free_vector(FF,     0, nXi);
  free_vector(DD,     0, nXi);
  free_vector(RR,     0, nXi);
  
  return 0;
}
