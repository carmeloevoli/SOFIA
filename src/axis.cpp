#include "axis.h"

size_t Axis::getLowerIndex(double x) const {
	assert(x >= front());
	size_t i = 0;
	while (get(i + 1) < x)
		i++;
	return i;
}

size_t Axis::getNearestIndex(double x) const {
	if (x < front())
		return 0;
	else if (x > back()) {
		return _N - 1;
	} else {
		size_t i_found = -1;
		double min_distance = 1e100;
		for (size_t i = 0; i < _N; ++i) {
			double distance = std::abs(x - get(i));
			if (distance < min_distance) {
				i_found = i;
				min_distance = distance;
			}
		}
		return i_found;
	}
}

void Axis::buildLinAxis(const double& min, const double& max, const size_t& size) {
	_N = size;
	const double dx = (max - min) / (double) (size - 1);
	for (size_t i = 0; i < size; ++i) {
		_axis.push_back(min + dx * i);
	}
}

void Axis::buildLogAxis(const double& min, const double& max, const size_t& size) {
	_N = size;
	const double delta_log = std::exp(std::log(max / min) / (size - 1));
	for (size_t i = 0; i < size; ++i) {
		_axis.push_back(std::exp(std::log(min) + (double) i * std::log(delta_log)));
	}
}
