#ifndef INCLUDE_SOURCES_H_
#define INCLUDE_SOURCES_H_

#include "axis.h"
#include "grid.h"
#include "params.h"
#include "pid.h"

using InjectionParams = std::pair<double, double>;

class Sources: public Grid {
	PID _pid;
	Axis _T;
	Axis _z;

public:
	Sources(const Axis& T, const Axis& z, const double& units) :
			Grid(T.size(), z.size(), units), _T(T), _z(z) {
	}

	const Axis& getAxisZ() const {
		return _z;
	}

	const Axis& getAxisT() const {
		return _T;
	}

	const PID& pid() const {
		return _pid;
	}

	void build(const PID& pid, const Params& params, const InjectionParams& injectionParams);

	double integrate() const;
};

#endif /* INCLUDE_SOURCES_H_ */
