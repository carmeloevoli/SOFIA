#ifndef INCLUDE_AXIS_H_
#define INCLUDE_AXIS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>
#include <vector>

class Axis {
protected:
	std::vector<double> _axis;
	size_t _N;
	double _units;

public:
	Axis() :
			_N(0), _units(0) {
	}

	Axis(double units) :
			_N(0), _units(units) {
	}

	virtual ~Axis() {
		clearAxis();
	}

	void resize(size_t N) {
		this->_N = N;
		_axis.resize(N);
	}

	/** Accessor / Mutator */
	double &get(const size_t& i) {
		return _axis[i];
	}

	/* Min / Max */
	double max() const {
		return *max_element(_axis.begin(), _axis.end());
	}

	double min() const {
		return *min_element(_axis.begin(), _axis.end());
	}

	/* Front / Back */
	double front() const {
		return _axis.front();
	}

	double back() const {
		return _axis.back();
	}

	double getUnits() const {
		return _units;
	}

	/** Accessor */
	const double &get(const size_t& i) const {
		return _axis[i];
	}

	double getValue(size_t i) const {
		return _axis[i];
	}

	/** Return a reference to the grid values */
	std::vector<double> &getAxis() {
		return _axis;
	}

	void clearAxis() {
		_axis.clear();
	}

	size_t size() const {
		return _N;
	}

	friend std::ostream& operator<<(std::ostream& stream, const Axis& axis) {
		stream << "(" << axis.min() / axis.getUnits();
		stream << "," << axis.max() / axis.getUnits();
		stream << "," << axis.size() << ")";
		return stream;
	}

	double getAverageStep() const {
		return std::abs(_axis.front() - _axis.back()) / (double) (_N - 1);
	}

	size_t getLowerIndex(double x) const;
	size_t getNearestIndex(double x) const;
	void buildLinAxis(const double& min, const double& max, const size_t& size);
	void buildLogAxis(const double& min, const double& max, const size_t& size);
};

#endif /* INCLUDE_AXIS_H_ */
