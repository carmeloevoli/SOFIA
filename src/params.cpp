#include "params.h"
#include <iostream>
#include <fstream>
#include <sstream>

Params::Params() {
	std::cout << "calling constructor of Params\n";
}

Params::~Params() {
	std::cout << "calling destructor of Params\n";
}

void Params::set_params(const std::string& key, const double& value) {
	if (key == "H")
		_H = value * cgs::kpc;
	else if (key == "dt")
		_dt = value * cgs::kyr;
	else if (key == "u")
		_u = value * cgs::km / cgs::sec;
	else if (key == "v_A")
		_v_A = value * cgs::km / cgs::sec;
	else if (key == "slope")
		_injection_slope = value;
	else if (key == "phi")
		_modulation_potential = value * cgs::GeV;
	else if (key == "id")
		_id = (int) value;
}

void Params::set_from_file(const std::string& filename) {
	std::ifstream infile(filename.c_str());
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string key;
		double value;
		if (!(iss >> key >> value)) {
			break;
		} // error
		set_params(key, value);
	}
}

void Params::print() {
	std::cout << "W_0     : " << _W_0 / (cgs::pc) << " pc\n";
	std::cout << "delta_k : " << _delta_k << "\n";
	std::cout << "k_0     : " << _k_0 * cgs::pc << " pc-1\n";
	std::cout << "B_0     : " << _B_0 / cgs::muG << " muG\n";
	std::cout << "u       : " << _u / (cgs::km / cgs::sec) << " km/s\n";
	std::cout << "slope   : " << _injection_slope << "\n";
	std::cout << "phi     : " << _modulation_potential / cgs::GeV << " GeV\n";
	std::cout << "H       : " << _H / cgs::kpc << " kpc\n";
	std::cout << "z_size  : " << _z_size << "\n";
	std::cout << "k_min   : " << _k_min / (1. / cgs::pc) << " pc-1\n";
	std::cout << "k_max   : " << _k_max / (1. / cgs::pc) << " pc-1\n";
	std::cout << "k_size  : " << _k_size << "\n";
	std::cout << "E_min   : " << _E_min / cgs::GeV << " GeV\n";
	std::cout << "E_max   : " << _E_max / cgs::GeV << " GeV\n";
	std::cout << "E_size  : " << _E_size << "\n";
	std::cout << "dt      : " << _dt / cgs::kyr << " kyr\n";
}
