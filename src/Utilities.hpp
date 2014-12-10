/***************************************************************
2012-09-18: Utility functions: error handling and 
memory allocation.
***************************************************************/

#ifndef UTILITIES_H   // include guard
#define UTILITIES_H 1

#include <string> 
#include <fstream>   // file I/O
#include <sstream>
#include <iostream>

// Error messaging:
void error (const std::string message);
int warning (const std::string message);
// Random real generator:
double random(long *idum);
double gasdev(long *idum);
// Importing & Exporting data:
void ImportVecs(double **matriz, long length, long nvecs, const char *filename);
void PrintVecs(double **table, long nrows, long ncols, std::ostream *output = &std::cout, int offset=0);
// Formatting:
std::string ZeroPad(int num, int max);


/*** TEMPLATES ***/

const int NR_END=1;
// NR functions for creating and destroying arrays and matrices:
template<typename type>
type *vector(long nl, long nh) {
  type *v;

  v=(type *) new type [(nh-nl+1+NR_END)];
  if (!v) error("allocation failure in vector()");
  return v-nl+NR_END;
}
template<typename type>
void free_vector(type *v, long nl, long nh) {
  delete [] (v+nl-NR_END);
}

template<typename type>
type **matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  type **m;

  m = (type **) new type* [nrow+NR_END];
  if (!m) error("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  m[nrl]= (type *) new type [nrow*ncol+NR_END];
  if (!m[nrl]) error("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  return m;
}
template <typename type>
void free_matrix(type **m, long nrl, long nrh, long ncl, long nch) {
  delete [] (m[nrl]+ncl-NR_END);
  delete [] (m+nrl-NR_END);
}

template <typename type>
type ***tensor3(long n1i, long n1f, long n2i, long n2f, long n3i, long n3f) {
  long i,j, N1=n1f-n1i+1, N2=n2f-n2i+1, N3=n3f-n3i+1;
  type ***t;

  /* allocate pointers to pointers to rows */
  t=(type ***) new type** [N1+NR_END];
  if (!t) error("allocation failure 1 in tensor3()");
  t += NR_END;
  t -= n1i;
  
  /* allocate pointers to rows and set pointers to them */
  t[n1i]=(type **) new type* [N1*N2+NR_END];
  if (!t[n1i]) error("allocation failure 2 in tensor3()");
  t[n1i] += NR_END;
  t[n1i] -= n2i;

  /* allocate rows and set pointers to them */
  t[n1i][n2i]=(type *) new type [N1*N2*N3+NR_END];
  if (!t[n1i][n2i]) error("allocation failure 3 in tensor3()");
  t[n1i][n2i] += NR_END;
  t[n1i][n2i] -= n3i;

  for(j=n2i+1;j<=n2f;j++) t[n1i][j]=t[n1i][j-1]+N3;
  for(i=n1i+1;i<=n1f;i++) {
    t[i]=t[i-1]+N2;
    t[i][n2i]=t[i-1][n2i]+N2*N3;
    for(j=n2i+1;j<=n2f;j++) t[i][j]=t[i][j-1]+N3;
  }
  return t;
}
template <typename type>
void free_tensor3(type ***t, long n1i, long n1f, long n2i, long n2f, long n3i, long n3f) {
  delete [] (t[n1i][n2i]+n3i-NR_END);
  delete [] (t[n1i]+n2i-NR_END);
  delete [] (t+n1i-NR_END);
}


// Import table (matrix[1..nr][1..nc]) from file:
template <typename type>
type **LoadTable(std::string filename, long *nr, long *nc, int offset=0) {
  using std::ifstream;
  using std::string;
  using std::istringstream;
  using std::ostringstream;
  long nrows=0, ncols=0, i, j;
  ifstream file;
  istringstream inputline; ostringstream outputline;
  string word, phrase;
  type **table;
  
  // Open file
  file.open(filename.c_str());
  if (!file.is_open()) error("LoadTable: cannot open file "+filename);
  
  // Count lines and columns:
  getline(file,phrase);
  outputline << phrase;
  inputline.str(outputline.str());
  while (inputline >> word) ncols++;
  while(!file.eof()) {getline(file,phrase); nrows++;}
  std::cout<<"LoadTable will allocate "<<nrows<<" lines and "<<ncols<<" columns for file "<<filename<<std::endl;

  // Loading values to table:
  file.clear();
  file.seekg(0);
  table=matrix<type>(offset,nrows+offset-1,offset,ncols+offset-1);
  for (i=offset; i<nrows+offset; i++)
    for (j=offset; j<ncols+offset; j++)
      if (!(file >> table[i][j])) error("LoadTable: more data expected in file "+filename);
  if(file >> word) error("LoadTable: data was ignored in "+filename);
  *nr=nrows; *nc=ncols;

  file.close();
  return table;
}


// Import list (vector[1..nitems]) from file:
template <typename type>
type *LoadList(std::string filename, long *nitems, int offset=0) {
  using std::ifstream;
  using std::string;
  long n=0, i;
  ifstream file;
  string item;
  type *list;
  
  // Open file
  file.open(filename.c_str());
  if (!file.is_open()) error("LoadList: cannot open file.");
  
  // Count entries:
  while(file >> item) n++;
  std::cout<<"LoadList will allocate "<<n<<" entries for file "<<filename<<std::endl;

  // Loading values to table:
  file.clear();
  file.seekg(0);
  list=vector<type>(offset,n+offset-1);
  for (i=offset; i<n+offset; i++)
    if (!(file >> list[i])) error("LoadList: more data expected in file "+filename);
  if(file >> item) error("LoadList: data was ignored in "+filename);
  *nitems = n;

  file.close();
  return list;
}

// Print table:
template <typename type>
void PrintTable(type **table, long nrows, long ncols, std::ostream *output = &std::cout, int offset=0) {
  long i, j;
  
  (*output).setf(std::ios_base::showpoint);
  (*output).precision(6);
  for (i=offset; i<nrows+offset; i++) {
    for (j=offset; j<ncols+offset; j++) {
      (*output).width(10); *output << table[i][j] << " ";
    }
    *output << std::endl;
  }      
}

#endif
