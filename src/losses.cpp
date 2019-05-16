#include "losses.h"
#include "cgs.h"
#include "utils.h"

#define mec2 (cgs::electron_mass_c2)
#define mpc2 (cgs::proton_mass_c2)
#define beta2 (beta * beta)
#define gamma2 (gamma * gamma)

void EnergyLosses::build(const PID& pid, const Params& params) {
	double factor = 3. * cgs::sigma_th * pow2(pid.get_Z()) * cgs::electron_mass_c2 / 4.;
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		auto T = _T.get(iT);
		auto beta = beta_func(T);
		auto gamma = gamma_func(T);
		auto Q_max = pid.get_A() * (T + 2. * mpc2);
		Q_max /= 1. + pow2(pid.get_A() * mpc2 + mec2) / (2. * mec2 * T * pid.get_A());
		auto B_H = std::log(2. * mec2 * Q_max * beta2 * gamma2 / pow2(cgs::Is_H)) - 2. * beta2;
		auto B_He = std::log(2. * mec2 * Q_max * beta2 * gamma2 / pow2(cgs::Is_He)) - 2. * beta2;
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			auto n_H = params.n_H0 * std::exp(-pow2(_z.get(iz) / params.h_d) / 2.);
			get(iT, iz) = factor / beta2 * n_H * (B_H + cgs::f_He * B_He);
		}
	}
}
