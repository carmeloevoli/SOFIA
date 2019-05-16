#include "advection.h"

void AdvectionVelocity::build_velocity_profile(const Params& params) {
	for (size_t i = 0; i < _z.size(); ++i) {
		auto z = _z.get(i);
		get(i) = params.u * std::tanh(z / params.h_d);
	}
}
