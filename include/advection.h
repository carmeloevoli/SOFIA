#ifndef INCLUDE_ADVECTION_H_
#define INCLUDE_ADVECTION_H_

#include <cmath>
#include "axis.h"
#include "cgs.h"
#include "grid.h"
#include "params.h"

class AdvectionVelocity: public Grid {
	Axis _z;

public:
	AdvectionVelocity(const Axis& z, const double& units) :
			Grid(z.size(), units), _z(z) {
	}

	void build_velocity_profile(const Params& params);
};

#endif /* INCLUDE_ADVECTION_H_ */
