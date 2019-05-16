#include "diffusion.h"
#include "utils.h"
#include <iostream>

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
