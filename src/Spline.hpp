#ifndef SPLINE_H    // include guard.
#define SPLINE_H 1

#include "FieldsDatabase.hpp" // Specifically for init using corrlnfields parameters.

class Spline {
private:
  int dim, N1, N2;
  double *x1, *x2, *f1d, **f2d, *d2fdx2, **d2fdxdy;
  void reset();
  void clear();
  void dealloc();
public:
  Spline();
  Spline(double *dom, double *f, int n);
  Spline(double *dom1, double *dom2, double **f, int n1, int n2);
  void init(double *x, double *y, int n);
  void init(double *dom1, double *dom2, double **f, int n1, int n2);
  void init(const FZdatabase & fieldlist, double ***Cov, int field, int l); // Specifically for init using corrlnfields parameters.
  double operator()(double x) const;
  double operator()(double x1, double x2) const;
  ~Spline();
};


// Standalone functions that serve as bakcbone (from Numerical Recipes)
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);


#endif
