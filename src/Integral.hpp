#ifndef INTEGRAL_H
#define INTEGRAL_H 1

#include "Cosmology.hpp"
typedef const Cosmology & param;

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double trapzd(double (*func)(param,double), double a, double b, int n, param p);
double qromb(double (*func)(param,double), double a, double b, param p);
double trapzd(double (*func)(double, double), double a, double b, int n, double p);
double qromb(double (*func)(double, double), double a, double b, double p);
double trapzd(double (*func)(double, double, param), double a, double b, int n, double z0, param p);
double qromb(double (*func)(double, double, param), double a, double b, double z0, param p);
double trapzd5(double (*func)(param, double, double), double a, double b, int n, double z0, param p);
double qromb5(double (*func)(param, double, double), double a, double b, double z0, param p);
#endif
