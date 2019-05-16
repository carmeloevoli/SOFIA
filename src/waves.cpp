#include "waves.h"
#include <iostream>

void WaveSpectrum::build_initial_condition(const Params& params) {
	for (size_t i = 0; i < _k.size(); ++i)
		for (size_t j = 0; j < _z.size(); ++j) {
			auto k = _k.get(i);
			auto x = k / params.k_0;
			get(i, j) = (x > 1) ? params.W_0 * std::pow(x, -params.delta_k) : 0.;
		}
	std::cout << " ... building power spectrum with <delta B> / B = " << params.W_0 * params.k_0 / (params.delta_k - 1.) << "\n";
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
