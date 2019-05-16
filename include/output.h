#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include <string>
#include "advection.h"
#include "axis.h"
#include "diffusion.h"
#include "losses.h"
#include "particle.h"
#include "sources.h"

class OutputManager {
	std::string _init;

public:
	OutputManager(std::string init_filename) :
			_init(init_filename) {
	}

	void dumpTimescales(const Sources& Q, const DiffusionCoefficient& D_zz, const MomentumDiffusionCoefficient & D_pp,
			const EnergyLosses& dpdt, const AdvectionVelocity& u) const;
	void dumpSpectrum(const Sources& Q, const DiffusionCoefficient& D, const EnergyLosses& dpdt) const;
	void dumpProfile(const Sources& Q, const DiffusionCoefficient& D, const EnergyLosses& dpdt,
			const AdvectionVelocity& u) const;
	void dumpSolution(const Particle& f) const;
};

#endif /* INCLUDE_OUTPUT_H_ */
