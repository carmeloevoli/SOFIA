#include "axis.h"

Axis::Axis() {
}

Axis::Axis(HaloSize H, GridSize N, GridStep dz) {
	grid.resize(N.get());
	int i = -1;
	std::generate(grid.begin(), grid.end(), [H, dz, i]() mutable {++i; return -H.get() + i * dz.get();});
	central_idx = (N.get() % 2 == 0) ? N.get() / 2 : (N.get() - 1) / 2;
}

double Axis::dz() const {
	return std::fabs(grid.at(1) - grid.at(0));
}

size_t Axis::get_central() const {
	return central_idx;
}

