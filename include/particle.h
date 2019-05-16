#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include "axis.h"
#include "grid.h"
#include "pid.h"

class Particle: public Grid {
	PID _pid;
	Axis _T;
	Axis _z;

public:
	Particle(const Axis& T, const Axis& z, const PID& pid, const double& units) :
			Grid(T.size(), z.size(), units), _T(T), _z(z), _pid(pid) {
		std::cout << "calling constructor of Particle\n";
	}

	~Particle() {
		std::cout << "calling destructor of Particle\n";
	}

	PID pid() const {
		return _pid;
	}

//	double get_z(const size_t& i) const {
//		return _z.get(i);
//	}
//
//	double get_T(const size_t& i) const {
//		return _T.get(i);
//	}
//
	const Axis& getAxisZ() const {
		return _z;
	}

	const Axis& getAxisT() const {
		return _T;
	}

	double integrate() const;
};

#endif /* INCLUDE_PARTICLE_H_ */
