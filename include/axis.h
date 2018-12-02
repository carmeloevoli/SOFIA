#ifndef AXIS_H_
#define AXIS_H_

#include <algorithm>
#include <cassert>
#include <cmath>

#include "common.h"
#include "grid.h"

class Axis: public Grid {
protected:
	size_t central_idx = 0;
public:
	Axis();
	Axis(HaloSize H, GridSize N, GridStep dz);
	double dz() const;
	size_t get_central() const;
};

#endif /* AXIS_H_ */
