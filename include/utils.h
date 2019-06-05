#ifndef INCLUDE_UTILS_H_
#define INCLUDE_UTILS_H_

#include <vector>
#include "axis.h"
#include "pid.h"

#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))

double beta_func(double T);
double gamma_func(double T);
double momentum_func(double T, size_t A);
double larmor_radius(double T, double B, int A, int Z);
double pres_func(double k, double B, int Z);
double Gamma_Integral(double slope);
double PowerlawSpectrum(double T_min, double T, double alpha);
double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);
double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);
double Gaussian(double x, double sigma);
Axis build_momentum_axis(const Axis& E, const PID& pid);

#endif /* INCLUDE_UTILS_H_ */
