#include "FieldsDatabase.hpp"
#include "FieldsDatabase.cpp"
#include <iostream>
#include "Utilities.hpp"
#include "Utilities.cpp"


void memtest(int *f, int *z, int Nf) {
  FZdatabase test(f,z,Nf);
}


int main () {
  int *f;
  int *z;
  int Nf, i, j, fname, zname, fi, zi;
  int *wrapper[2];
  long nr, nc;
  std::string filename;
  
  filename.assign("test.dat");
  LoadVecs(wrapper, filename, &nr, &nc, 0, 1);
  printf("vai passar para vecs individuais...\n");
  f = wrapper[0];
  z = wrapper[1];
  Nf = (int)nr;

  // Test for memory leackage:
  printf("mem test...\n");
  for (i=0; i<1000000; i++) memtest(f, z, Nf);
  printf("done.\n");
  
  FZdatabase fzinfo(f, z, Nf);

  printf("\nNfs: %d   Nzs: %d   Nfields: %d\n\n", fzinfo.Nfs(), fzinfo.Nzs(), fzinfo.Nfields());

  printf("LOOP field (redshift): \n");
  for(i=0; i<fzinfo.Nfs(); i++) {
    for(j=0; j<fzinfo.Nz4f(i); j++) {
      fzinfo.fFixedName(i, j, &fname, &zname);
      printf("%d  %d, ", fname, zname);
    }
    printf("\n");
  }

  printf("\nLOOP redshift (field): \n");
  for(i=0; i<fzinfo.Nzs(); i++) {
    for(j=0; j<fzinfo.Nf4z(i); j++) {
      fzinfo.zFixedName(j, i, &fname, &zname);
      printf("%d  %d, ", fname, zname);
    }
    printf("\n");
  }

  printf("Index for fFixed:\n");
  for(i=0; i<fzinfo.Nfields(); i++) {
    fzinfo.Index2fFixed(i, &fi, &zi);
    printf("%d  %d\n", fi, zi);
  }
  printf("Index for zFixed:\n");
  for(i=0; i<fzinfo.Nfields(); i++) {
    fzinfo.Index2zFixed(i, &fi, &zi);
    printf("%d  %d\n", fi, zi);
  }
  
  return 0;
}
