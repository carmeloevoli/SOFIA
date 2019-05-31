#include "diffusion.h"
#include "utils.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#define LIMIT 1000
#define EPSREL 1e-5

void DiffusionCoefficient::build_QLT(const PID& pid, const WaveSpectrum& W, const double& B_0) {
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		auto T = _T.get(iT);
		auto v = beta_func(T) * cgs::c_light;
		auto rL = larmor_radius(T, B_0, pid.get_A(), pid.get_Z());
		auto D_B = v * rL / 3.0;
		auto k = 1. / rL;
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			double kW = k * W.getInterpolated(k, iz);
			auto value = D_B / kW;
			get(iT, iz) = value;
		}
	}
}

double compute_integral_qag(gsl_integration_workspace * w, gsl_function * F, double x_lo, double x_hi) {
	double result, error;
	int key = 3;
	gsl_integration_qag(F, x_lo, x_hi, 0, EPSREL, LIMIT, key, w, &result, &error);
	return result;
}

struct gslDiffusionClassParams {
	double rL;
	size_t iz;
	const WaveSpectrum* W;
};

double gslDiffusionClassFunction(double logkres, void * p) {
	double kres = std::exp(logkres);
	gslDiffusionClassParams params = *(gslDiffusionClassParams *) p;
	double rL = params.rL;
	double rL_kres = rL * kres;
	double kres_W = kres * params.W->getInterpolated(kres, params.iz);
	return 1. / rL_kres * (1. - 1. / pow2(rL_kres)) / kres_W;
}

void DiffusionCoefficient::build_Skilling(const PID& pid, const Params& params, const WaveSpectrum& W) {
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		auto T = _T.get(iT);
		auto v = beta_func(T) * cgs::c_light;
		auto rL = larmor_radius(T, params.B_0, pid.get_A(), pid.get_Z());
		auto D_B = v * rL / pow2(2. * M_PI);
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			gsl_function F;
			gslDiffusionClassParams gslParams = { rL, iz, &W };
			F.params = &gslParams;
			F.function = &gslDiffusionClassFunction;
			gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
			assert(params.k_max > 1. / rL);
			double I = compute_integral_qag(w, &F, std::log(1. / rL), std::log(params.k_max));
			gsl_integration_workspace_free(w);
			auto value = D_B * I;
			get(iT, iz) = value;
		}
	}
}

void MomentumDiffusionCoefficient::build(const PID& pid, const Params& params, const DiffusionCoefficient& D_zz) {
	double delta = 2. - params.delta_k; // TODO Problem with D_xx generic!
	double w = 1.; // TODO compute from W?
	double factor = 4. * pow2(params.v_A);
	factor /= 3 * delta * (4 - pow2(delta)) * (4 - delta) * w;
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		double p = momentum_func(_T.get(iT), pid.get_A());
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			double value = factor * p * p;
			get(iT, iz) = value / D_zz.get(iT, iz);
		}
	}
}

double gslMomentumDiffusionClassFunction(double logkres, void * p) {
	double kres = std::exp(logkres);
	gslDiffusionClassParams params = *(gslDiffusionClassParams *) p;
	double rL = params.rL;
	double rL_kres = rL * kres;
	double kres_W = kres * params.W->getInterpolated(kres, params.iz);
	return 1. / rL_kres * (1. - 1. / pow2(rL_kres)) * kres_W;
}

void MomentumDiffusionCoefficient::build_Skilling(const PID& pid, const Params& params, const WaveSpectrum& W) {
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		auto T = _T.get(iT);
		auto v = beta_func(T) * cgs::c_light;
		auto rL = larmor_radius(T, params.B_0, pid.get_A(), pid.get_Z());
		auto p = momentum_func(T, pid.get_A());
		auto factor = pow2(M_PI * params.v_A * p) / (v * rL);
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			gsl_function F;
			gslDiffusionClassParams gslParams = { rL, iz, &W };
			F.params = &gslParams;
			F.function = &gslMomentumDiffusionClassFunction;
			gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
			assert(params.k_max > 1. / rL);
			double I = compute_integral_qag(w, &F, std::log(1. / rL), std::log(params.k_max));
			gsl_integration_workspace_free(w);
			auto value = factor * I;
			get(iT, iz) = value;
		}
	}
}
