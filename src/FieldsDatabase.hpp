#ifndef FZDATABASE_H    // include guard.
#define FZDATABASE_H 1

#include <vector>
#include <cstddef> // For NULL pointer.
#include <string>

struct entry {
  int IDnum;
  int Ncopies;
  std::vector<int> index;
};

class FZdatabase {
private:
  // Field names storage:
  int *fullF, *fullZ;                        // For an index 'n' store the fields and redshift bin names here.
  std::vector<entry>f, z;                    // Individual lists for every field or z, used for looping.
  int **ifsubz;
  int **izsubf;
  int Nf, Nz, Nfield;                       // Number of fields (small f), number of redshift bins, total number of Fields (Nf*Nz)
  // Fields info:
  double *means, *shifts, **zrange;
  int *ftypes;
  // Auxiliary functions:
  int GetPos(int ent, std::vector<entry> & list);
  void RegisterEntry(int *full, int i, std::vector<entry> & vec, int **ixsuby);
  void Init(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0);
public:
  // Structural functions:
  FZdatabase();
  FZdatabase(const std::string & filename);
  FZdatabase(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0);
  void Load(const std::string & filename);
  void Build(int *fullF0, int *fullZ0, int Nfield0, int *ftype0, double **zrange0, double *mean0, double *shift0);
  ~FZdatabase();
  // Data size functions:
  int Nzs() const;                                // Return the number of distinct redshift bins, no matter the field (small f).
  int Nfs() const;                                // Return the number of fields (small f).
  int Nfields() const;                            // Return total number of Fields (Nf*Nz).
  int Nz4f(int fi) const;                         // Return number of redshift bins for a particular field index fi.
  int Nf4z(int zi) const;                         // Return number of fields found at a particular redshift bin index zi.
  // Looping functions:
  int  zFixedIndex(int fi, int zi, int *n=NULL) const; // To loop inside a fixed redshift.
  int  fFixedIndex(int fi, int zi, int *n=NULL) const; // To loop inside a fixed field.
  // Field name functions:
  void Index2fFixed(int n, int *fi, int *zi) const;
  void Index2zFixed(int n, int *fi, int *zi) const;
  int  Name2Index(int fName, int zName, int *n=NULL, bool warn=1) const;
  void Index2Name(int n,          int *fName, int *zName) const;
  void zFixedName(int fi, int zi, int *fName, int *zName) const;
  void fFixedName(int fi, int zi, int *fName, int *zName) const;
  // Field info functions:
  double mean(int n)  const;
  double shift(int n) const;
  int    ftype(int n) const;
  double zmin(int n)  const;
  double zmax(int n)  const;
  // Field organization information:
  int CheckZ4Int() const;
  // Extract two field and redshift NAME pairs from string:
  void String2NamePair(const std::string & str, int *f1, int *z1, int *f2, int *z2) const;
};


#endif
