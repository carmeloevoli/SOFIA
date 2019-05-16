#ifndef INCLUDE_WAVES_H_
#define INCLUDE_WAVES_H_

#include <cmath>
#include "axis.h"
#include "grid.h"
#include "params.h"

class WaveSpectrum: public Grid {
	Axis _k;
	Axis _z;

public:
	WaveSpectrum(const Axis& k, const Axis& z, const double& units) :
			Grid(k.size(), z.size(), units), _k(k), _z(z) {
	}

	void build_initial_condition(const Params& params);
	double getInterpolated(const double& k, const size_t iz) const;
};

#endif /* INCLUDE_WAVES_H_ */
