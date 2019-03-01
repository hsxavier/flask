#include "Utilities.hpp"
#include "FieldsDatabase.hpp"
#include "definitions.hpp" // For ftype possibilities.
#include <string>

/************************************/
/*** Internal auxiliary functions ***/
/************************************/


// Return the position of entry in list or -1 of not found:
int FZdatabase::GetPos(int ent, std::vector<entry> & list) {
  int i;
  for (i=0; i<list.size(); i++) if (ent==list[i].IDnum) return i;
  return -1;
}


// Either create new entry in 'vec' or register it as another appearance of previous one:
void FZdatabase::RegisterEntry(int *full, int i, std::vector<entry> & vec, int **ixsuby) {
  int pos;
  entry aux;

  pos = GetPos(full[i], vec);
  // If entry is new, add to list of distinct values:
  if (pos==-1) {
    aux.IDnum    = full[i];
    aux.Ncopies  = 1;
    aux.index.resize(1); aux.index[0]=i;
    vec.push_back(aux);
    ixsuby[i][0] = vec.size()-1;
    ixsuby[i][1] = 0;
  }
  // If not, add the index to the entries list:
  else { 
    vec[pos].index.push_back(i); 
    vec[pos].Ncopies++; 
    ixsuby[i][0] = pos;
    ixsuby[i][1] = vec[pos].Ncopies-1;
  }
}


// Initialize an empty FZdatabase:
void FZdatabase::Init(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0) {
  int i;

  // Allocate memory for field and redshift list from input:
  Nfield = Nfield0;
  fullF  = vector<int>(0, Nfield-1);
  fullZ  = vector<int>(0, Nfield-1);
  ifsubz = matrix<int>(0, Nfield-1, 0, 1);
  izsubf = matrix<int>(0, Nfield-1, 0, 1);
  // Field info:
  ftypes = vector<int>   (0, Nfield-1);
  zrange = matrix<double>(0, Nfield-1, 0, 1);
  means  = vector<double>(0, Nfield-1);
  shifts = vector<double>(0, Nfield-1);


  // Copy list to internal memory, organizing in database:
  for (i=0; i<Nfield; i++) {
    fullF[i] = fullF0[i]; 
    fullZ[i] = fullZ0[i];
    RegisterEntry(fullF, i, f, ifsubz);
    RegisterEntry(fullZ, i, z, izsubf);
    ftypes[i]    = ftype0[i];
    zrange[i][0] = zrange0[i][0];
    zrange[i][1] = zrange0[i][1];
    means[i]     = mean0[i];
    shifts[i]    = shift0[i];
  }
  Nf = f.size();
  Nz = z.size();  
}


/************************/
/*** Public functions ***/
/************************/


// Empty constructor:
FZdatabase::FZdatabase() {
  Nf     = 0;
  Nz     = 0;
  Nfield = 0;  
}


// Constructor with input:
FZdatabase::FZdatabase(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0) {
  Init(fullF0, fullZ0, Nfield0, ftype0, zrange0, mean0, shift0);
}
// Constructor with input:
FZdatabase::FZdatabase(const std::string & filename) {
  Nf     = 0;
  Nz     = 0;
  Nfield = 0;   
  Load(filename);
}


void FZdatabase::Load(const std::string & filename) {
  using namespace definitions;
  char message[200];
  int i, zi, *fName, *zName;
  long long1, long2;
  double **aux;
  
  // Safe guard against double loading:
  if (Nfield>0) error("FZdatabase.Load: second time loading not implemented.");

  // Load info file:
  aux    = LoadTable<double>(filename, &long1, &long2);
  // Field names:
  Nfield = (int)long1;
  fullF  = vector<int>   (0, Nfield-1);
  fullZ  = vector<int>   (0, Nfield-1);
  ifsubz = matrix<int>   (0, Nfield-1, 0, 1);
  izsubf = matrix<int>   (0, Nfield-1, 0, 1);
  // Field info:
  means  = vector<double>(0, Nfield-1);
  shifts = vector<double>(0, Nfield-1);
  ftypes = vector<int>   (0, Nfield-1);
  zrange = matrix<double>(0, Nfield-1, 0, 1);
  
  // Parse information to separate arrays:
  for (i=0; i<Nfield; i++) {
    // Field names:
    fullF[i]     = (int)aux[i][0]; 
    fullZ[i]     = (int)aux[i][1];
    RegisterEntry(fullF, i, f, ifsubz);
    RegisterEntry(fullZ, i, z, izsubf);
    // Field info:
    means[i]     =      aux[i][2];  
    shifts[i]    =      aux[i][3];    
    ftypes[i]    = (int)aux[i][4];
    zrange[i][0] =      aux[i][5];
    zrange[i][1] =      aux[i][6]; 
  }
  free_matrix(aux, 0, long1-1, 0, long2-1);
  Nf = f.size();
  Nz = z.size(); 
 
  // A few checks on the input:
  for (i=0; i<Nfield; i++) {
    if (zrange[i][0]>zrange[i][1])    warning("FZdatabase.Load: zmin > zmax for a field.");
    if (ftypes[i]!=fgalaxies && ftypes[i]!=flensing) warning("FZdatabase.Load: unknown field type in FIELDS_INFO file.");
    if(means[i]+shifts[i]<=0) {
      sprintf(message, "FZdatabase.Load: mean+shift at position %d must be greater than zero (LOGNORMAL realizations only).", i); 
      warning(message);
    }
  }
  // Check if all redshift bins with the same name have the same redshift range:
  for (zi=0; zi<Nz; zi++) 
    for (i=1; i<z[zi].Ncopies; i++) 
      if (zrange[z[zi].index[i]][0] != zrange[z[zi].index[0]][0] || zrange[z[zi].index[i]][1] != zrange[z[zi].index[0]][1])
	warning("FZdatabase.Load: redshift bins with same index should have same zmin and zmax, but they do not.");
}

// Build fields and redshifts database erasing previous one:
void FZdatabase::Build(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0) {
  // Destroy previous allocation if existent:
  if (Nf>0) { free_vector(fullF, 0, Nf); f.resize(0); }
  if (Nz>0) { free_vector(fullZ, 0, Nz); z.resize(0); }
  if (Nfield>0) {
    free_vector(ftypes, 0, Nfield-1);
    free_matrix(zrange, 0, Nfield-1, 0, 1);
    free_vector(means,  0, Nfield-1);
    free_vector(shifts, 0, Nfield-1);
  }
  // Initialize with new information:
  Init(fullF0, fullZ0, Nfield0, ftype0, zrange0, mean0, shift0);
}


// Destructor:
FZdatabase::~FZdatabase() {
  if (Nf>0)     free_vector(fullF,  0, Nf);
  if (Nz>0)     free_vector(fullZ,  0, Nz);
  if (Nfield>0) free_vector(means,  0, Nfield-1);
  if (Nfield>0) free_vector(shifts, 0, Nfield-1);
  if (Nfield>0) free_vector(ftypes, 0, Nfield-1);
  if (Nfield>0) free_matrix(zrange, 0, Nfield-1, 0, 1);
  if (Nfield>0) free_matrix(ifsubz, 0, Nfield-1, 0, 1);
  if (Nfield>0) free_matrix(izsubf, 0, Nfield-1, 0, 1);
}


// Return number of field types:
int FZdatabase::Nfs() const {
  return Nf;
}


// Return number of redshift bins:
int FZdatabase::Nzs() const {
  return Nz;
}


// Return total number of Fields (Nf*Nz):
int FZdatabase::Nfields() const {
  return Nfield;
}


// Return number of redshift bins for a particular field index 'fi':
int FZdatabase::Nz4f(int fi) const {
  if (fi>Nf || fi<0) error("FZdatabase.Nz4f: requested field does not exist.");
  return f[fi].Ncopies;
}


// Return number of fields found at a particular redshift bin index 'zi':
int FZdatabase::Nf4z(int zi) const {
  if (zi>Nz || zi<0) error("FZdatabase.Nf4z: requested redshift bin does not exist.");
  return z[zi].Ncopies; 
}


// Given a Field index 'n', returns the field index 'fi' and the redshift sub index 'zi':
void FZdatabase::Index2fFixed(int n, int *fi, int *zi) const {
  *fi = ifsubz[n][0];
  *zi = ifsubz[n][1];
}


// Given a Field index 'n', returns the redshift index 'zi' and the field sub index 'fi':
void FZdatabase::Index2zFixed(int n, int *fi, int *zi) const {
  *fi = izsubf[n][1];
  *zi = izsubf[n][0];
}
  


// Return the index of the Field with z index 'zi' and f sub-index 'fi':
// (This is for LOOPS over f with z fixed)
int FZdatabase::zFixedIndex(int fi, int zi, int *n) const {
  if (zi<0 || zi>Nz)            error("FZdatabase.zFixedIndex: unknown redshift index.");
  if (fi<0 || fi>z[zi].Ncopies) error("FZdatabase.zFixedIndex: unknown field sub-index.");
  if (n!=NULL) *n = z[zi].index[fi]; 
  return z[zi].index[fi];
}


// Return the index of the Field with f index 'fi' and z sub-index 'zi':
// (This is for LOOPS over z with f fixed)
int FZdatabase::fFixedIndex(int fi, int zi, int *n) const {
  if (fi<0 || fi>Nf)            error("FZdatabase.fFixedIndex: unknown field index.");
  if (zi<0 || zi>f[fi].Ncopies) error("FZdatabase.fFixedIndex: unknown redshift sub-index.");
  if (n!=NULL) *n = f[fi].index[zi]; 
  return f[fi].index[zi];
}


// Return field NAME and redshift bin NAME for an index n:
void FZdatabase::Index2Name(int n, int *fName, int *zName) const {
  *fName = fullF[n];
  *zName = fullZ[n];
}


// Return index of a Field given a name: 
int FZdatabase::Name2Index(int fName, int zName, int *n, bool warn) const {
  const int failed=-1;
  int i;

  for (i=0; i<Nfield; i++) 
    if (fName==fullF[i] && zName==fullZ[i]) {
      if (n!=NULL) *n = i;
      return i;
    }
  
  if (warn) warning("FZdatabase.Name2Index: could not find requested field.");  
  *n = failed;
  return failed;
}


// Return field NAME and redshift bin NAME for an z index 'zi' and f sub-index 'fi':
void FZdatabase::zFixedName(int fi, int zi, int *fName, int *zName) const {
  Index2Name(zFixedIndex(fi, zi), fName, zName);
}


// Return field NAME and redshift bin NAME for an f index 'fi' and z sub-index 'zi':
void FZdatabase::fFixedName(int fi, int zi, int *fName, int *zName) const {
  Index2Name(fFixedIndex(fi, zi), fName, zName);
}


// Return info for a Field index:
double FZdatabase::mean(int n)  const { return means[n];    }
double FZdatabase::shift(int n) const { return shifts[n];   }
int    FZdatabase::ftype(int n) const { return ftypes[n];   }
double FZdatabase::zmin(int n)  const { return zrange[n][0];}
double FZdatabase::zmax(int n)  const { return zrange[n][1];}


// Error checking for LoS integration (density fields must have continuous redshift coverage):
int FZdatabase::CheckZ4Int() const {
  using namespace definitions;
  int k=0, f, z, i, j;

  for (f=0; f<Nf; f++) {
    i = fFixedIndex(f, 0);
    if (ftypes[i]==fgalaxies) {
      k++; // Count density fields.
      for (z=1; z<Nz4f(f); z++) {
	fFixedIndex(f, z-1, &i); fFixedIndex(f, z, &j); 
	if (zmax(i) != zmin(j)) // Check if z bins are sequential and contiguous.
	  warning("FZdatabase.CheckZ4Int: expecting sequential AND contiguous redshift slices for galaxies");
      }
    }
  }
  return k;
}


// Extract 2 sets of field and redshift NAMES from a string in the format <A>f<F1>z<Z1>f<F2>z<Z2><B>:
void FZdatabase::String2NamePair(const std::string & str, int *f1, int *z1, int *f2, int *z2) const {
  int f1pos=-1, z1pos, f2pos, z2pos, end;
  bool StillLooking=true, IsNum1, IsNum2, IsNum3, IsNum4;

  // Find f<X>z<Y>f<W>z<Z> pattern:
  do {
    // Look for field and redshift name markers:
    f1pos = str.find('f', f1pos+1);
    z1pos = str.find('z', f1pos+1);
    f2pos = str.find('f', z1pos+1);
    z2pos = str.find('z', f2pos+1);
    if (f1pos==-1 || z1pos==-1 || f2pos==-1 || z2pos==-1) 
      error("FZdatabase.String2NamePair: could not find all field and redshift name markers \'f\' and \'z\'.");
    // Check if they follow the correct pattern (they are followed by numbers):
    IsNum1 = IsNumber(str.substr(f1pos+1,z1pos-f1pos-1));
    IsNum2 = IsNumber(str.substr(z1pos+1,f2pos-z1pos-1));
    IsNum3 = IsNumber(str.substr(f2pos+1,z2pos-f2pos-1));
    IsNum4 = IsNumber(str.substr(z2pos+1,1));
    if (IsNum1 && IsNum2 && IsNum3 && IsNum4) StillLooking=false; 
  } while (StillLooking);

  // Get field and redshift names:
  *f1 = str2int(str.substr(f1pos+1,z1pos-f1pos-1));
  *z1 = str2int(str.substr(z1pos+1,f2pos-z1pos-1));
  *f2 = str2int(str.substr(f2pos+1,z2pos-f2pos-1));
  for (end=z2pos+2; end<str.length() && isdigit(str[end])!=0; end++);
  *z2 = str2int(str.substr(z2pos+1,end-z2pos-1));
} 

