#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include "axis.h"
#include "grid.h"
#include "pid.h"

class Particle: public Grid {
	Axis _E;
	Axis _z;
	PID _pid;

public:
	Particle(const Axis& E, const Axis& z, const PID& pid, const double& units) :
			Grid(E.size(), z.size(), units), _E(E), _z(z), _pid(pid) {
		std::cout << "calling constructor of Particle\n";
	}

	~Particle() {
		std::cout << "calling destructor of Particle\n";
	}

	PID pid() const {
		return _pid;
	}

	double get_z(const size_t& i) const {
		return _z.get(i);
	}

	double get_E(const size_t& i) const {
		return _E.get(i);
	}

	const Axis& getAxisZ() const {
		return _z;
	}

	const Axis& getAxisE() const {
		return _E;
	}

	double integrate() const;
};

#endif /* INCLUDE_PARTICLE_H_ */
