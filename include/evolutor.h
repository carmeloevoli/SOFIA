#ifndef INCLUDE_EVOLUTOR_H_
#define INCLUDE_EVOLUTOR_H_

#include "advection.h"
#include "diffusion.h"
#include "losses.h"
#include "particle.h"
#include "sources.h"

//class TridiagMatrix {
//protected:
//	std::vector<double> _rhs;
//	std::vector<double> _central_diagonal;
//	std::vector<double> _upper_diagonal;
//	std::vector<double> _lower_diagonal;
//	std::vector<double> _fcr_up;
//	size_t _N;
//
//public:
//	TridiagMatrix() :
//			_N(0) {
//	}
//
//	TridiagMatrix(const size_t N) :
//			_N(N) {
//		_rhs.resize(N - 2);
//		_central_diagonal.resize(N - 2);
//		_upper_diagonal.resize(N - 3);
//		_lower_diagonal.resize(N - 3);
//		_fcr_up.resize(N - 2);
//	}
//
//	double& central_diagonal(size_t i) {
//		return _central_diagonal.at(i);
//	}
//
//	double& upper_diagonal(size_t i) {
//		return _upper_diagonal.at(i);
//	}
//
//	double& lower_diagonal(size_t i) {
//		return _lower_diagonal.at(i);
//	}
//
//	double& rhs(size_t i) {
//		return _rhs.at(i);
//	}
//
//	void fill(size_t i, double Li, double Ci, double Ui, double dt, double fi, double fup, double fdo, double Qi) {
//		double dt_half = 0.5 * dt;
//		_central_diagonal.at(i - 1) = 1. - dt_half * Ci;
//		if (i != _N - 2) {
//			_upper_diagonal.at(i - 1) = -dt_half * Ui;
//		}
//		if (i != 1) {
//			_lower_diagonal.at(i - 2) = -dt_half * Li;
//		}
//		_rhs.at(i - 1) = fi * (2. - _central_diagonal.at(i - 1));
//		_rhs.at(i - 1) += dt * Qi;
//		if (i != _N - 2) {
//			_rhs.at(i - 1) -= fup * _upper_diagonal.at(i - 1);
//		}
//		if (i != 1) {
//			_rhs.at(i - 1) -= fdo * _lower_diagonal.at(i - 2);
//		}
//	}
//
//};

class TimeEvolutor {
	double _dt = 0;
	size_t _max_counter = 0;
	bool evolveInSpace = false;
	bool evolveInMomentum = false;
	Axis _T;
	Axis _z;
	Axis _k;

public:
	TimeEvolutor(const Axis& T, const Axis& z, double dt, size_t max_counter) :
			_T(T), _z(z), _dt(dt), _max_counter(max_counter) {
		std::cout << " ... building time evolutor with max time : " << dt * max_counter / cgs::Myr << " Myr\n";
	}

	void doEvolveInSpace() {
		evolveInSpace = true;
	}

	void doEvolveInMomentum() {
		evolveInMomentum = true;
	}

	size_t num_operators() {
		size_t n = 0;
		n += (evolveInSpace) ? 1 : 0;
		n += (evolveInMomentum) ? 1 : 0;
		return n;
	}

	Particle evolve(const Sources& Q, const DiffusionCoefficient& D_zz, const MomentumDiffusionCoefficient& D_pp,
			const EnergyLosses& dpdt, const AdvectionVelocity& u_z, PID& pid);
	void evolve_in_z(const Sources& Q, const DiffusionCoefficient& D_zz, const AdvectionVelocity& u_z, Particle& f,
			const size_t& number_of_operators);
	void evolve_in_p(const Axis& p, const Sources& Q, const MomentumDiffusionCoefficient& D_pp,
			const EnergyLosses& dpdt, const AdvectionVelocity& u_z, Particle& f, const size_t& number_of_operators);

	void dumpSpectrum(const Particle& f, size_t idisc, size_t counter) const;
};

#endif /* INCLUDE_EVOLUTOR_H_ */
