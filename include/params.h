#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

#include <cmath>
#include <string>
#include "cgs.h"

class Params {
private:
	double _E_min = 1e-2 * cgs::GeV;
	double _E_max = 1e3 * cgs::GeV;
	size_t _E_size = 5 * 32;
	double _k_min = 1. / (1e-2 * cgs::pc);
	double _k_max = 1. / (1e-8 * cgs::pc);
	size_t _k_size = 6 * 16;
	double _H = 4. * cgs::kpc;
	size_t _z_size = 101;
	double _dt = cgs::kyr;
	double _u = 10. * cgs::km / cgs::sec;
	double _injection_slope = 4.4;
	double _modulation_potential = 0.7 * cgs::GeV;
	double _E_SN = 1e51 * cgs::erg;
	double _R_SN = 1. / (40. * cgs::year);
	double _R_d = 10. * cgs::kpc;
	double _h_d = 100. * cgs::pc;
	double _W_0 = 10. * cgs::pc;
	double _delta_k = 5. / 3.;
	double _k_0 = 1. / (100. * cgs::pc);
	double _B_0 = 1. * cgs::muG;
	double _n_H0 = 1. / cgs::cm_3;
	double _v_A = 30. * cgs::km / cgs::sec;
	size_t _id = 0;

public:
	Params();
	virtual ~Params();
	void print();
	void set_from_file(const std::string& filename);
	void set_params(const std::string& key, const double& value);

	const double& T_min = _E_min;
	const double& T_max = _E_max;
	const double& k_min = _k_min;
	const double& k_max = _k_max;
	const double& dt = _dt;
	const double& H = _H;
	const double& u = _u;
	const double& injection_slope = _injection_slope;
	const double& modulation_potential = _modulation_potential;
	const double& E_SN = _E_SN;
	const double& R_SN = _R_SN;
	const double& R_d = _R_d;
	const double& h_d = _h_d;
	const double& W_0 = _W_0;
	const double& k_0 = _k_0;
	const double& delta_k = _delta_k;
	const double& B_0 = _B_0;
	const double& n_H0 = _n_H0;
	const double& v_A = _v_A;
	const size_t& T_size = _E_size;
	const size_t& k_size = _k_size;
	const size_t& z_size = _z_size;
	const size_t& id = _id;
};

#endif /* INCLUDE_PARAMS_H_ */
