#include "cgs.h"
#include <cmath>
#include <gsl/gsl_integration.h>
#include "../include/utils.h"

#define LIMIT 1000

double beta_func(double T) {
	auto value = std::sqrt(T * (T + 2. * cgs::proton_mass_c2));
	value /= T + cgs::proton_mass_c2;
	return value;
}

double gamma_func(double T) {
	auto value = T / cgs::proton_mass_c2 + 1;
	return value;
}

double momentum_func(double T, size_t A) {
	auto value = beta_func(T) / cgs::c_light * (double) A * (T + cgs::proton_mass_c2);
	return value;
}

double larmor_radius(double T, double B, int A, int Z) {
	auto R = beta_func(T) * (T + cgs::proton_mass_c2) * ((double) A / (double) Z);
	return R / cgs::elementary_charge / B;
}

double gsl_Gamma_Integrand(double x, void * params) {
	double alpha = *(double *) params;
	double f = std::pow(x, 2. - alpha);
	f *= std::sqrt(x * x + 1) - 1;
	return f;
}

double Gamma_Integral(double slope) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
	double result, error;
	gsl_function F;
	F.function = &gsl_Gamma_Integrand;
	F.params = &slope;
	gsl_integration_qagiu(&F, 0, 0, 1e-6, LIMIT, w, &result, &error);
	gsl_integration_workspace_free(w);
	return 4. * M_PI * result;
}

double Gaussian(double x, double sigma) {
	return pow(2. * M_PI * pow2(sigma), -0.5) * std::exp(-pow2(x / sigma) / 2.);
}

double PowerlawSpectrum(double T_min, double T, double alpha) {
	double value = 0;
	if (T > T_min) {
		auto x = std::sqrt(pow2(T / cgs::proton_mass_c2 + 1) - 1);
		value = std::pow(x, -alpha);
	}
	return value;
}

Axis build_momentum_axis(const Axis& E, const PID& pid) {
	Axis p(cgs::GeV_c);
	p.resize(E.size());
	for (size_t i = 0; i < E.size(); ++i) {
		double E_total = pid.get_A() * (E.get(i) + cgs::proton_mass_c2);
		p.get(i) = beta_func(E.get(i)) * E_total / cgs::c_light;
	}
	return p;
}

double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y, const double& x_new) {
	auto value = double();
	if (x_new > x.front() && x_new <= x.back()) {
		size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
		double t = (x_new - x.at(i - 1)) / (x.at(i) - x.at(i - 1));
		value = y.at(i - 1) * (1. - t) + y.at(i) * t;
	}
	return value;
}

double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y, const double& x_new) {
	auto value = double();
	if (x_new >= x.front() && x_new <= x.back()) {
		size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
		double t = std::log(x_new) - std::log(x.at(i - 1));
		t /= std::log(x.at(i)) - std::log(x.at(i - 1));
		value = std::log(y.at(i - 1)) * (1. - t) + std::log(y.at(i)) * t;
		value = std::exp(value);
	}
	return value;
}

#undef LIMIT
