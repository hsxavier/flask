#ifndef FZDATABASE_H    // include guard.
#define FZDATABASE_H 1

#include <vector>
#include <cstddef> // For NULL pointer.


struct entry {
  int IDnum;
  int Ncopies;
  std::vector<int> index;
};

class FZdatabase {
private:
  // Information storage:
  int *fullF, *fullZ;                       // For an index 'n' store the fields and redshift bin names here.
  std::vector<entry>f, z;                   // Individual lists for every field or z, used for looping.
  int Nf, Nz, Nfield;                       // Number of fields (small f), number of redshift bins, total number of Fields (Nf*Nz)
  // Auxiliary functions:
  int GetPos(int ent, std::vector<entry> & list);
  void RegisterEntry(int *full, int i, std::vector<entry> & vec);
  void Init(int *fullF0, int *fullZ0, int Nfield0);
public:
  // Structural functions:
  FZdatabase();
  FZdatabase(int *fullF0, int *fullZ0, int Nfield0);
  void Build(int *fullF0, int *fullZ0, int Nfield0);
  ~FZdatabase();
  // Data size functions:
  int Nzs();                                // Return the number of distinct redshift bins, no matter the field (small f).
  int Nfs();                                // Return the number of fields (small f).
  int Nfields();                            // Return total number of Fields (Nf*Nz).
  int Nz4f(int fi);                         // Return number of redshift bins for a particular field index fi.
  int Nf4z(int zi);                         // Return number of fields found at a particular redshift bin index zi.
  // Looping functions:
  int  zFixedIndex(int fi, int zi, int *n = NULL); // To loop inside a fixed redshift.
  int  fFixedIndex(int fi, int zi, int *n = NULL); // To loop inside a fixed field.
  // Field name functions:
  int  Name2Index(int fName, int zName, int *n = NULL);
  void Index2Name(int n,          int *fName, int *zName);
  void zFixedName(int fi, int zi, int *fName, int *zName);
  void fFixedName(int fi, int zi, int *fName, int *zName);
};


#endif
