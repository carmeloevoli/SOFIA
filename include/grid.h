#ifndef INCLUDE_GRID_H_
#define INCLUDE_GRID_H_

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

class Grid {
protected:
	std::vector<double> _grid;
	size_t _Nx, _Nz;
	double _units;

public:
	Grid() :
			_Nx(0), _Nz(0), _units(0) {
	}

	Grid(size_t Nx, size_t Nz, double units) :
			_units(units) {
		resize(Nx, Nz);
	}

	Grid(size_t Nz, double units) :
			_units(units) {
		resize(1, Nz);
	}

	virtual ~Grid() {
		clearGrid();
	}

	void resize(size_t NE, size_t Nz) {
		this->_Nx = NE;
		this->_Nz = Nz;
		_grid.resize(NE * Nz);
	}

	/** Accessor / Mutator */
	double &get(size_t ix, size_t iz) {
		return _grid[ix * _Nz + iz];
	}

	double &get(const size_t& i) {
		return _grid[i];
	}

	/* Min / Max */
	double max() const {
		return *max_element(_grid.begin(), _grid.end());
	}

	double min() const {
		return *min_element(_grid.begin(), _grid.end());
	}

	double getUnits() const {
		return _units;
	}

	/** Accessor */
	const double &get(size_t iE, size_t iz) const {
		return _grid[iE * _Nz + iz];
	}

	const double &get(const size_t& i) const {
		return _grid[i];
	}

	double getValue(size_t iE, size_t iz) const {
		return _grid[iE * _Nz + iz];
	}

	/** Return a reference to the grid values */
	std::vector<double> &getGrid() {
		return _grid;
	}

	void clearGrid() {
		_grid.clear();
	}

	size_t getNx() const {
		return _Nx;
	}

	size_t getNz() const {
		return _Nz;
	}

	size_t size() const {
		return (_Nx * _Nz);
	}

	friend std::ostream& operator<<(std::ostream& stream, const Grid& grid) {
		stream << "(" << grid.min() / grid.getUnits();
		stream << "," << grid.max() / grid.getUnits();
		stream << "," << grid.size() << ")";
		return stream;
	}
};

#endif /* INCLUDE_GRID_H_ */
