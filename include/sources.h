#ifndef INCLUDE_SOURCES_H_
#define INCLUDE_SOURCES_H_

#include "axis.h"
#include "grid.h"
#include "params.h"
#include "pid.h"

using InjectionParams = std::pair<double, double>;

class Sources: public Grid {
	Axis _E;
	Axis _z;
	PID _pid;

public:
	Sources(const Axis& E, const Axis& z, const double& units) :
			Grid(E.size(), z.size(), units), _E(E), _z(z) {
	}

	const Axis& getAxisZ() const {
		return _z;
	}

	const Axis& getAxisT() const {
		return _E;
	}

	const PID& pid() const {
		return _pid;
	}

	void build(const PID& pid, const Params& params, const InjectionParams& injectionParams);

	double integrate() const;
};

#endif /* INCLUDE_SOURCES_H_ */
