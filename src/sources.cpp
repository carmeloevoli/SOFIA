#include "sources.h"
#include "utils.h"

void Sources::build(const PID& pid, const Params& params, const InjectionParams& injectionParams) {
	_pid = pid;
	double Mc = pid.get_A() * cgs::proton_mass_c;
	double disk_surface = 4. * M_PI * pow2(params.R_d);
	double Q_0 = injectionParams.first * params.E_SN * params.R_SN;
	Q_0 /= disk_surface * cgs::c_light * pow4(Mc) * Gamma_Integral(injectionParams.second);
	for (size_t iT = 0; iT < _T.size(); ++iT) {
		auto T = _T.get(iT);
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			auto z = _z.get(iz);
			get(iT, iz) = Q_0 * Gaussian(z, params.h_d) * PowerlawSpectrum(0.1 * cgs::GeV, T, injectionParams.second); // TODO min rigidity
		}
	}
}

double Sources::integrate() const {
	double value = 0;
	auto dz = _z.getAverageStep();
	auto p = build_momentum_axis(_T, _pid);
	for (size_t iz = 0; iz < _z.size(); ++iz) {
		double I_p = 0;
		for (size_t iT = 0; iT < _T.size() - 1; ++iT) {
			auto fup = pow2(p.get(iT + 1)) * _T.get(iT + 1) * get(iT + 1, iz);
			auto fdo = pow2(p.get(iT)) * _T.get(iT) * get(iT, iz);
			I_p += (p.get(iT + 1) - p.get(iT)) * 0.5 * (fup + fdo);
		}
		value += I_p;
	}
	return 4. * M_PI * value * dz;
}
