#include "sources.h"
#include "utils.h"

void Sources::build(const PID& pid, const Params& params, const InjectionParams& injectionParams) {
	_pid = pid;
	double Mc = pid.get_A() * cgs::proton_mass_c;
	double disk_surface = 4. * M_PI * pow2(params.R_d);
	double Q_0 = injectionParams.first * params.E_SN * params.R_SN;
	Q_0 /= disk_surface * cgs::c_light * pow4(Mc) * Gamma_Integral(injectionParams.second);
	for (size_t iE = 0; iE < _E.size(); ++iE) {
		auto E = _E.get(iE);
		for (size_t iz = 0; iz < _z.size(); ++iz) {
			auto z = _z.get(iz);
			get(iE, iz) = Q_0 * Gaussian(z, params.h_d) * PowerlawSpectrum(0.1 * cgs::GeV, E, injectionParams.second); // TODO min rigidity
		}
	}
}

double Sources::integrate() const {
	double value = 0;
	auto dz = _z.getAverageStep();
	auto p = build_momentum_vector(_E, _pid);
	for (size_t iz = 0; iz < _z.size(); ++iz) {
		double I_p = 0;
		for (size_t iE = 0; iE < _E.size() - 1; ++iE) {
			auto fup = pow2(p.get(iE + 1)) * _E.get(iE + 1) * get(iE + 1, iz);
			auto fdo = pow2(p.get(iE)) * _E.get(iE) * get(iE, iz);
			I_p += (p.get(iE + 1) - p.get(iE)) * 0.5 * (fup + fdo);
		}
		value += I_p;
	}
	return 4. * M_PI * value * dz;
}
