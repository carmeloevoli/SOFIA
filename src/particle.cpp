#include "particle.h"
#include "utils.h"

double Particle::integrate() const {
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
