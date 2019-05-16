#include "evolutor.h"
#include "tridiag.h"
#include "utils.h"

void TimeEvolutor::dumpSpectrum(const Particle& f, size_t idisc, size_t counter) const {
	std::string filename = "output/f_spectrum_at_" + std::to_string(counter) + ".txt";
	std::cout << "printing f on this file: " << filename << " ... ";
	std::ofstream outfile(filename.c_str());
	outfile << "# T [GeV] \t f [] \n";
	outfile << std::scientific;
	for (size_t i = 0; i < _T.size(); ++i) {
		outfile << _T.get(i) / _T.getUnits() << "\t";
		outfile << f.get(i, idisc) / f.getUnits() << "\t";
		outfile << "\n";
	}
	outfile.close();
	std::cout << "... done!" << "\n";
}

Particle TimeEvolutor::evolve(const Sources& Q, const DiffusionCoefficient& D_zz,
		const MomentumDiffusionCoefficient& D_pp, const EnergyLosses& dpdt, const AdvectionVelocity& u_z, PID& pid) {
	auto p = build_momentum_axis(_T, pid);
	size_t counter = 0;
	size_t idisc = _z.getNearestIndex(0.);
	Particle f(_T, _z, pid, 1. / cgs::cm_3 / cgs::GeV_c_3);
	auto I_Q = Q.integrate();
	while (counter < _max_counter) {
		if (evolveInSpace)
			evolve_in_z(Q, D_zz, u_z, f, num_operators());
		if (evolveInMomentum)
			evolve_in_p(p, Q, D_pp, dpdt, u_z, f, num_operators());
		if (counter % 1000 == 0) {
			auto t = _dt * (double) counter;
			auto I_f = f.integrate();
			std::cout << t / cgs::Myr << " " << I_Q * t / cgs::erg << " " << I_f / cgs::erg << " " << " "
					<< I_f / I_Q / cgs::Myr << "\n";
		}
		//dumpSpectrum(f, idisc, counter);
		counter++;
	}
	return f;
}

void TimeEvolutor::evolve_in_p(const Axis& p, const Sources& Q, const MomentumDiffusionCoefficient& D_pp,
		const EnergyLosses& dpdt, const AdvectionVelocity& u, Particle& f, const size_t& number_of_operators) {
	double dt_half = 0.5 * _dt;
	double dz = _z.getAverageStep();
	size_t p_size = p.size(), z_size = _z.size();
	std::vector<double> rhs(p_size - 2);
	std::vector<double> central_diagonal(p_size - 2);
	std::vector<double> upper_diagonal(p_size - 3);
	std::vector<double> lower_diagonal(p_size - 3);
	std::vector<double> fcr_up(p_size - 2);

	for (size_t iz = 1; iz < z_size - 1; ++iz) {
		double dudz = (u.get(iz + 1) - u.get(iz - 1)) / 2. / dz;
		for (size_t ip = 1; ip < p_size - 1; ++ip) {
			double b_i = dudz / 3.0 / (p.get(ip + 1) / p.get(ip) - 1.0);
			double Lp = 0, Cp = 0, Up = 0;
			Cp += b_i;
			Up += b_i;

			Cp += dpdt.get(ip, iz) / p.get(ip) * (1. / (p.get(ip + 1) / p.get(ip) - 1.0) - 2.);
			Up += dpdt.get(ip + 1, iz) / (p.get(ip + 1) - p.get(ip));

			double p_c = p.get(ip + 1) - p.get(ip - 1);
			double p_up = p.get(ip + 1) - p.get(ip);
			double p_do = p.get(ip) - p.get(ip - 1);
			double Dpp_c = D_pp.get(ip, iz);
			double Dpp_up = D_pp.get(ip + 1, iz);
			double Dpp_do = D_pp.get(ip - 1, iz);
			double p_i = p.get(ip);

			Lp += 2. * Dpp_c / pow2(p_c) * (p_c / p_up + 0.5 * (Dpp_up - Dpp_do) / Dpp_c + p_c / p_i);
			Cp += 2. * Dpp_c / pow2(p_c) * (p_c / p_up + p_c / p_do);
			Up += 2. * Dpp_c / pow2(p_c) * (p_c / p_do - 0.5 * (Dpp_up - Dpp_do) / Dpp_c - p_c / p_i);

			central_diagonal.at(ip - 1) = 1. + dt_half * Cp;
			if (ip != 1) {
				lower_diagonal.at(ip - 2) = -dt_half * Lp;
			}
			if (ip != p_size - 2) {
				upper_diagonal.at(ip - 1) = -dt_half * Up;
			}
			rhs.at(ip - 1) = f.get(ip, iz) * (2. - central_diagonal.at(ip - 1));
			rhs.at(ip - 1) += _dt * Q.get(ip, iz) / (double) number_of_operators;
			if (ip != p_size - 2) {
				rhs.at(ip - 1) -= f.get(ip + 1, iz) * upper_diagonal.at(ip - 1);
			}
			if (ip != 1) {
				rhs.at(ip - 1) -= f.get(ip - 1, iz) * lower_diagonal.at(ip - 2);
			}
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (size_t ip = 1; ip < p_size - 1; ++ip) {
			double value = fcr_up.at(ip - 1);
			f.get(ip, iz) = std::max(value, 0.);
		}
		f.get(0, iz) = f.get(1, iz);
	} // for
}

void TimeEvolutor::evolve_in_z(const Sources& Q, const DiffusionCoefficient& D_zz, const AdvectionVelocity& u_z,
		Particle& f, const size_t& number_of_operators) {
	double dt_half = 0.5 * _dt;
	double dz = _z.getAverageStep();
	size_t T_size = _T.size(), z_size = _z.size();
	size_t idisc = _z.getNearestIndex(0.);
	std::vector<double> rhs(z_size - 2);
	std::vector<double> central_diagonal(z_size - 2);
	std::vector<double> upper_diagonal(z_size - 3);
	std::vector<double> lower_diagonal(z_size - 3);
	std::vector<double> fcr_up(z_size - 2);
	for (size_t iT = 0; iT < T_size - 1; iT++) {
		for (size_t iz = 1; iz < z_size - 1; iz++) {
			double UZ = 0, CZ = 0, LZ = 0;
			const double Dz = D_zz.get(iT, iz);
			const double DzUp = (iz < z_size - 1) ? D_zz.get(iT, iz + 1) : Dz;
			const double DzDo = (iz > 0) ? D_zz.get(iT, iz - 1) : Dz; // TODO check this
			const double Dz_dz2 = Dz / pow2(dz);
			const double dDz_4_dz2 = (DzUp - DzDo) / 4. / pow2(dz);

			UZ += Dz_dz2 + dDz_4_dz2;
			CZ += 2. * Dz_dz2;
			LZ += Dz_dz2 - dDz_4_dz2;

			if (iz > idisc) {
				LZ += u_z.get(iz - 1) / dz;
				CZ += u_z.get(iz) / dz;
				UZ += 0;
			} else if (iz == idisc) {
				LZ += 0;
				CZ += 0;
				UZ += 0;
			} else {
				LZ += 0;
				CZ += -u_z.get(iz) / dz;
				UZ += -u_z.get(iz + 1) / dz;
			}

			central_diagonal.at(iz - 1) = 1. + dt_half * CZ;
			if (iz != z_size - 2) {
				upper_diagonal.at(iz - 1) = -dt_half * UZ;
			}
			if (iz != 1) {
				lower_diagonal.at(iz - 2) = -dt_half * LZ;
			}

			rhs.at(iz - 1) = f.get(iT, iz) * (2. - central_diagonal.at(iz - 1));
			rhs.at(iz - 1) += _dt * Q.get(iT, iz) / (double) number_of_operators;
			if (iz != z_size - 2) {
				rhs.at(iz - 1) -= f.get(iT, iz + 1) * upper_diagonal.at(iz - 1);
			}
			if (iz != 1) {
				rhs.at(iz - 1) -= f.get(iT, iz - 1) * lower_diagonal.at(iz - 2);
			}
		}

		gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, fcr_up);

		for (size_t iz = 1; iz < z_size - 1; ++iz) {
			double value = 0.5 * (fcr_up.at(iz - 1) + fcr_up.at(z_size - 2 - iz));
			f.get(iT, iz) = (value > 0) ? value : 0.;
		}
	} // for iE
} // evolve_in_z
