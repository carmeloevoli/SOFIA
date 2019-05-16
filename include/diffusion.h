#ifndef INCLUDE_DIFFUSION_H_
#define INCLUDE_DIFFUSION_H_

#include "axis.h"
#include "grid.h"
#include "pid.h"
#include "waves.h"

class DiffusionCoefficient: public Grid {
	Axis _T;
	Axis _z;

public:
	DiffusionCoefficient(const Axis& T, const Axis& z, const double& units) :
			Grid(T.size(), z.size(), units), _T(T), _z(z) {
	}

	void build_QLT(const PID& pid, const WaveSpectrum& W, const double& B_0);
};

class MomentumDiffusionCoefficient: public Grid {
	Axis _T;
	Axis _z;

public:
	MomentumDiffusionCoefficient(const Axis& T, const Axis& z, const double& units) :
			Grid(T.size(), z.size(), units), _T(T), _z(z) {
	}

	void build(const PID& pid, const Params& params, const DiffusionCoefficient& D);
};

#endif /* INCLUDE_DIFFUSION_H_ */
