#include <iostream>
#include <cstdlib>              // For function 'system'.
#include "Utilities.hpp"        // Error handling, tensor allocations.
#include "interpol.h"           // Interpolation.

#define  MAXINPUT 500
#define  NCLMAX   500


/********************/
/*** Main Program ***/
/********************/
int main (int argc, char *argv[]) {
  using std::cout; using std::endl; using std::string; using std::ofstream; // Basic stuff.
  using std::cin;
  char mes1[MAXINPUT], mes2[MAXINPUT];
  std::string filename;
  std::ofstream outfile;                                                    // File for output.
  std::ifstream infile, covfile;
  int i, j, n, a1, a2, b1, b2, N1, N2, Nentries[NCLMAX], llout=5;
  bool **IsSet;
  long ncols, Nl, maxNl;
  double ***ll, ***Cov, *wrapper[2], **CovMatrix;
  void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2);
  void CountEntries(std::string filename, long *nr, long *nc);

  // Listing files to use, given in command line.
  mes1[0]='\0';
  mes2[0]='\0';
  for(i=1; i<argc; i++) sprintf(mes1,"%s %s", mes1, argv[i]);
  sprintf(mes2, "ls%s > gencovl.temp", mes1);
  system(mes2);

  // Get file list and find out how many C(l)s there are:
  i=0; N1=0; N2=0, maxNl=0;
  infile.open("gencovl.temp");
  if (!infile.is_open()) error("Cannot open file gencovl.temp");
  while (infile >> filename) {
    getcovid(filename, &a1, &a2, &b1, &b2);
    if (a1>N1) N1=a1; if (b1>N1) N1=b1;     // Get number of fields.
    if (a2>N2) N2=a2; if (b2>N2) N2=b2;     // Get number of z bins.
    CountEntries(filename, &Nl, &ncols);    // Get number of Nls.
    if (ncols!=2) error("Wrong number of columns in file "+filename);
    if (Nl>maxNl) maxNl=Nl;                 // Record maximum number of ls.
    Nentries[i] = Nl;                       // Record number of ls for each file.
    i++;
    if (i>NCLMAX) error("Reached maximum number of C(l)s. Increase NCLMAX.");
  }
  infile.clear();
  infile.seekg(0);
  cout << "Nfields: " << N1 << " Nzs: " << N2 << endl;
  
  // Allocate memory to store C(l)s:
  // First two indexes are CovMatrix indexes and last is for ll.
  ll  = tensor3<double>(1, N1*N2, 1, N1*N2, 0, maxNl);
  Cov = tensor3<double>(1, N1*N2, 1, N1*N2, 0, maxNl);
  IsSet = matrix<bool>(1, N1*N2, 1, N1*N2); // For bookkepping.
  for(i=1; i<=N1*N2; i++) for(j=1; j<=N1*N2; j++) IsSet[i][j]=0;
    
  // Read C(l)s and store in data-cube:
  n=0;
  while (infile >> filename) {
    // Find CovMatrix indexes of C(l):
    getcovid(filename, &a1, &a2, &b1, &b2);
    i=(a1-1)*N1+a2; j=(b1-1)*N1+b2;
    cout << filename << " goes to ["<<i<<", "<<j<<"]" << endl;
    // Import data:
    wrapper[0] = &(ll[i][j][0]);
    wrapper[1] = &(Cov[i][j][0]);
    ImportVecs(wrapper, Nentries[n], 2, filename.c_str());
    IsSet[i][j]=1;
  };
  infile.close();
  
  // Set the rest of CovMatrix as a symmetric matrix:
  for(i=1; i<=N1*N2; i++) 
    for(j=1; j<=N1*N2; j++)
      if (IsSet[i][j]==0) {
	for (n=0; n<maxNl; n++) {
	  ll[i][j][n]=ll[j][i][n];
	  Cov[i][j][n]=Cov[j][i][n];
	}
	IsSet[i][j]=IsSet[j][i];
	if (IsSet[i][j]==0) error ("Cannot set full CovMatrix, missing C(l).");
      }
  
  // Interpolate and save covariance matrices for each ll:
  CovMatrix = matrix<double>(1, N1*N2, 1, N1*N2);
  for (n=1; n<=llout; n++) {
    sprintf(mes1, "cov-ll%d.dat",n);
    outfile.open(mes1);
    for (i=1; i<=N1*N2; i++)
      for (j=1; j<=N1*N2; j++)
	CovMatrix[i][j] = Interpol(&(ll[i][j][0]), maxNl, &(Cov[i][j][0]), (double)n);
    PrintTable(CovMatrix, N1*N2, N1*N2, &outfile, 1);
    outfile.close();
  }
    
  // End of program:
  system("rm -f gencovl.temp");
  cout << "\nTotal number of warnings: " << warning("count") << endl;
  cout<<endl;
  return 0;
}


// Get four numbers separated by characters that specify the fields and redshifts
// of the correlation function.
void getcovid(const std::string filename, int *a1, int *a2, int *b1, int *b2) {
  int i=0, num, index, fileL;
  
  fileL=filename.length();
  // LOOP over the four indexes that indentifies the C(l):
  for(index=1; index<=4; index++) {
    num=0;
    // Find a number:
    while (isdigit(filename.c_str()[i])==0) {i++; if(i>=fileL && index!=4) error("getcovid: cannot find four numbers.");}
    // Read the number:
    while (isdigit(filename.c_str()[i])!=0) {num = num*10 + (filename.c_str()[i]-'0'); i++;}
    // Save it as an index:
    switch (index) {
    case 1: *a1 = num; break;
    case 2: *a2 = num; break;
    case 3: *b1 = num; break;
    case 4: *b2 = num; break;
    }
  }
  // Check if there are more numbers in filename:
  while (i<=fileL) {
    if (isdigit(filename.c_str()[i])!=0) error("getcovid: found more numbers than expected.");
    i++;
  }
}


// Find out number of columns and rows in file:
void CountEntries(std::string filename, long *nr, long *nc) {
  using std::ifstream;
  using std::string;
  using std::istringstream;
  using std::ostringstream;
  long nrows=0, ncols=0;
  ifstream file;
  istringstream inputline; ostringstream outputline;
  string word, phrase;
  
  // Open file
  file.open(filename.c_str());
  if (!file.is_open()) error("CountEntries: cannot open file.");
  
  // Count lines and columns:
  getline(file,phrase);
  outputline << phrase;
  inputline.str(outputline.str());
  while (inputline >> word) ncols++;
  while(!file.eof()) {getline(file,phrase); nrows++;}

  file.close();
  *nr=nrows+1;
  *nc=ncols;
}
