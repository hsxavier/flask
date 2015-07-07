#include <cmath>


// Compute the lognormal distribution shift parameter from the distribution's mean 'm', variance 'v' and skewness 'g':
double Moments2Shift(double m, double v, double g) {
  const double OneThird = 1.0/3.0; 
  double g2, y;

  g2 = g*g;
  y  = pow( (2.0 + g2 + g*sqrt(4.0+g2)) / 2.0, OneThird);

  return sqrt(v)*(1.0 + 1.0/y + y)/g - m;
} 


// Returns the mean of the normal distribution associated with a lognormal one with mean 'm', variance 'v' and shift:
double gmu(double m, double v, double shift) {
  return log((m + shift)/sqrt(1.0 + v/pow(m + shift, 2.0)));
}


// Returns the mean of the normal distribution associated with a lognormal one with mean 'm', variance 'v' and shift:
double gsigma(double m, double v, double shift) {
  return sqrt(log(1.0 + v/pow(m + shift, 2.0)));
}
