#ifndef INCLUDE_PID_H_
#define INCLUDE_PID_H_

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

class PID {
public:
	PID() {
		set(0, 0);
	}

	PID(const int& Z, const int& A) {
		assert(A > 0);
		assert(Z <= A);
		set(Z, A);
	}

	virtual ~PID() {
	}

	void set(const int& Z, const int& A) {
		_Z = Z;
		_A = A;
		_id = A * 1000 + Z;
	}

	int get_Z() const {
		return _Z;
	}

	int get_A() const {
		return _A;
	}

	double get_Z_over_A() const {
		return (_A > 0) ? fabs((double) _Z / (double) _A) : 0;
	}

	double get_A_over_Z() const {
		return (_A > 0) ? fabs((double) _A / (double) _Z) : 0;
	}

	int get_id() const {
		return _id;
	}

	bool operator==(const PID &other) const {
		return _id == other._id;
	}

	bool operator!=(const PID &other) const {
		return _id != other._id;
	}

	bool operator<(const PID &other) const {
		return _id < other._id;
	}

	bool operator>(const PID &other) const {
		return _id > other._id;
	}

	friend std::ostream& operator<<(std::ostream& stream, const PID& pid) {
		stream << "(" << pid.get_A() << "," << pid.get_Z() << ")";
		return stream;
	}

	std::string to_string() const {
		std::string ss;
		ss = "(" + std::to_string(_Z) + "," + std::to_string(_A) + ")";
		return ss;
	}

protected:
	int _Z;
	int _A;
	int _id;
};

typedef std::pair<PID, PID> Channel;

static const PID H1 = PID(1, 1);
static const PID H2 = PID(1, 2);
static const PID He3 = PID(2, 3);
static const PID He4 = PID(2, 4);

#endif /* INCLUDE_PID_H_ */
