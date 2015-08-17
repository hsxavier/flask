#include "FieldsDatabase.hpp"
#include "FieldsDatabase.cpp"
#include <iostream>
#include "Utilities.hpp"
#include "Utilities.cpp"

int main () {
  int *f;
  int *z;
  int Nf, i, j, fname, zname;
  int **wrapper;
  long nr, nc;
  std::string filename;
  
  filename.assign("test.dat");
  LoadVecs(wrapper, filename, &nr, &nc, 0, 1);
  f = wrapper[0];
  z = wrapper[1];
  Nf = (int)nr;

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


  return 0;
}
