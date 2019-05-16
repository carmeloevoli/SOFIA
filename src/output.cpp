#include "output.h"
#include "utils.h"
#include <fstream>

void OutputManager::dumpTimescales(const Sources& Q, const DiffusionCoefficient& D_zz,
		const MomentumDiffusionCoefficient& D_pp, const EnergyLosses& dpdt, const AdvectionVelocity& u) const {
	auto z = Q.getAxisZ();
	auto T = Q.getAxisT();
	auto pid = Q.pid();
	double H = z.max();
	size_t idisc = z.getNearestIndex(0.);
	std::string filename = _init + "_timescales.txt";
	std::ofstream fout(filename.c_str());
	fout << std::scientific;
	for (size_t i = 0; i < T.size(); ++i) {
		fout << T.get(i) / T.getUnits() << " ";
		fout << beta_func(T.get(i)) << " ";
		double p = beta_func(T.get(i)) / cgs::c_light * pid.get_A() * (T.get(i) + cgs::proton_mass_c2);
		fout << H * H / 2. / D_zz.get(i, idisc) / cgs::Myr << " ";
		fout << p * p / D_pp.get(i, idisc) / cgs::Myr << " ";
		fout << H / u.max() / cgs::Myr << " ";
		fout << p / dpdt.get(i, idisc) / cgs::Myr << " ";
		fout << "\n";
	}
	fout.close();
}

void OutputManager::dumpSpectrum(const Sources& Q, const DiffusionCoefficient& D, const EnergyLosses& dpdt) const {
	auto z = Q.getAxisZ();
	auto T = Q.getAxisT();
	size_t idisc = z.getNearestIndex(0.);
	std::string filename = _init + "_galaxy_spectrum.txt";
	std::ofstream fout(filename.c_str());
	fout << std::scientific;
	for (size_t iT = 0; iT < T.size(); ++iT) {
		fout << T.get(iT) / T.getUnits() << " ";
		fout << Q.get(iT, idisc) / Q.getUnits() << " ";
		fout << D.get(iT, idisc) / D.getUnits() << " ";
		fout << dpdt.get(iT, idisc) / dpdt.getUnits() << " ";
		fout << "\n";
	}
	fout.close();
}

void OutputManager::dumpProfile(const Sources& Q, const DiffusionCoefficient& D, const EnergyLosses& dpdt,
		const AdvectionVelocity& u) const {
	auto z = Q.getAxisZ();
	auto T = Q.getAxisT();
	size_t iT = T.getNearestIndex(10. * cgs::GeV);
	std::string filename = _init + "_galaxy_profile.txt";
	std::ofstream fout(filename.c_str());
	fout << std::scientific;
	for (size_t iz = 0; iz < z.size(); ++iz) {
		fout << z.get(iz) / z.getUnits() << " ";
		fout << Q.get(iT, iz) / Q.getUnits() << " ";
		fout << D.get(iT, iz) / D.getUnits() << " ";
		fout << dpdt.get(iT, iz) / dpdt.getUnits() << " ";
		fout << u.get(iz) / u.getUnits() << " ";
		fout << "\n";
	}
	fout.close();
}

void OutputManager::dumpSolution(const Particle& f) const {
	auto z = f.getAxisZ();
	auto T = f.getAxisT();
	auto A = f.pid().get_A();
	size_t idisc = z.getNearestIndex(0.);
	std::string filename = _init + "_solution.txt";
	std::ofstream fout(filename.c_str());
	fout << std::scientific;
	for (size_t iT = 0; iT < T.size(); ++iT) {
		double T_n = T.get(iT);
		double E_total = A * (T_n + cgs::proton_mass_c2);
		double p = beta_func(T_n) * E_total / cgs::c_light;
		fout << T_n / cgs::GeV << " ";
		fout << A * p * p * f.get(iT, idisc) / (1. / cgs::GeV / cgs::m_2 / cgs::sec) << " ";
		fout << "\n";
	}
	fout.close();
}

