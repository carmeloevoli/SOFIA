#ifndef INCLUDE_LOSSES_H_
#define INCLUDE_LOSSES_H_

#include "axis.h"
#include "grid.h"
#include "params.h"
#include "pid.h"

class EnergyLosses: public Grid {
	Axis _T;
	Axis _z;

public:
	EnergyLosses(const Axis& T, const Axis& z, const double& units) :
			Grid(T.size(), z.size(), units), _T(T), _z(z) {
	}

	void build(const PID& pid, const Params& params);
};

#endif /* INCLUDE_LOSSES_H_ */
