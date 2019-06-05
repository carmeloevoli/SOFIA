#include "waves.h"
#include "utils.h"
#include <iostream>
#include <gsl/gsl_integration.h>

void WaveSpectrum::build_initial_condition(const Params& params) {
	for (size_t i = 0; i < _k.size(); ++i)
		for (size_t j = 0; j < _z.size(); ++j) {
			auto k = _k.get(i);
			auto x = k / params.k_0;
			get(i, j) = (x > 1) ? params.W_0 * std::pow(x, -params.delta_k) : 0.;
		}
	auto deltab_b = params.W_0 * params.k_0 / (params.delta_k - 1.);
	std::cout << " ... building power spectrum with <delta B> / B = " << deltab_b << "\n";
}

double WaveSpectrum::getInterpolated(const double& k, const size_t iz) const {
	if (k < _k.front() || k > _k.back())
		return 0;
	else {
		size_t const i = _k.getLowerIndex(k);
		double t = std::log(k) - std::log(_k.get(i));
		t /= std::log(_k.get(i + 1)) - std::log(_k.get(i));
		double v = std::log(get(i, iz)) * (1. - t) + std::log(get(i + 1, iz)) * t;
		return std::exp(v);
	}
}

double f_H(double p) {
	double f_0 = 9.2e-10 / pow3(cgs::GeV_c * cgs::meter);
	return f_0 * std::pow(p / (10. * cgs::GeV_c), -4.8);
}

double cyclotron_damping_integral_func(double log_p, void * params) {
	double p = std::exp(log_p);
	return p * p * f_H(p);
}

double cyclotron_damping_integral(double p_min, double p_max) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

	double result, error;
	double alpha = 1.0;

	gsl_function F;
	F.function = &cyclotron_damping_integral_func;
	F.params = &alpha;

	gsl_integration_qag(&F, std::log(p_min), std::log(p_max), 0, 1e-7, 3, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

double cyclotron_damping_integral_Brunetti_func(double log_p, void * params) {
	double k = *(double *) params;
	double p = std::exp(log_p);
	const double h = cgs::MeV_c;
	double dfdp = -f_H(p + 2 * h) + 8. * f_H(p + h) - 8. * f_H(p - h) + f_H(p - 2. * h);
	dfdp /= 12. * h;
	double r_L = p * cgs::c_light / cgs::elementary_charge / (1. * cgs::muG);
	double mu_alpha = 1. / k / r_L;
	return -p * p * p * dfdp * (1. - pow2(mu_alpha));
}

double cyclotron_damping_integral_Brunetti(double p_min, double p_max, double k) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

	double result, error;

	gsl_function F;
	F.function = &cyclotron_damping_integral_Brunetti_func;
	F.params = &k;

	gsl_integration_qag(&F, std::log(p_min), std::log(p_max), 0, 1e-7, 3, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

void WaveSpectrum::sandbox(const Params& params) const {
//	for (double p = cgs::GeV_c; p < cgs::TeV_c; p *= 1.1) {
//		double E = std::sqrt(pow2(p * cgs::c_light) + pow2(cgs::proton_mass_c2));
//		double beta = p * cgs::c_light / E;
//		double I = beta * cgs::c_light * p * p * f_H(p);
//		std::cout << std::scientific << p / cgs::GeV_c << "\t" << I / (1. / cgs::GeV_c / cgs::m_2 / cgs::sec) << "\n";
//	}
//	exit(1);

	for (size_t i = 0; i < _k.size(); ++i) {
		double k = _k.get(i);
		double Z = 1; // TODO generalize to other than protons
		double v_A = params.v_A;
		double v_A_over_c = v_A / cgs::c_light;
		double p_res = pres_func(k, params.B_0, Z);
		double p_max = cgs::PeV_c; // TODO put infinity?
		double I = cyclotron_damping_integral(p_res, p_max);
		double cyclotron_damping = 2 * pow2(M_PI * Z * cgs::elementary_charge * v_A_over_c) / k * I;

		double I_B = cyclotron_damping_integral_Brunetti(p_res, p_max, k);

		double cyclotron_damping_Brunetti = pow2(2. * M_PI * Z * cgs::elementary_charge * v_A_over_c) / k * I_B;

		double t_c = 1. / v_A / k / std::sqrt(k * get(i, 10));

		std::cout << std::scientific << p_res / cgs::GeV_c << "\t" << k * cgs::pc << "\t" << 1. / cyclotron_damping / cgs::Myr
				<< "\t" << 1. / cyclotron_damping_Brunetti / cgs::Myr << "\t" << t_c / cgs::Myr << "\n";
	}
}
