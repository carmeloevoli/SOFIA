#include "particle.h"
#include "utils.h"
#include "cgs.h"

double Particle::integrate() const {
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
