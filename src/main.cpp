#include <iomanip>
#include <iostream>
#include <vector>

#include "axis.h"
#include "grid.h"
#include "params.h"
#include "particle.h"
#include "advection.h"
#include "diffusion.h"
#include "evolutor.h"
#include "losses.h"
#include "output.h"
#include "sources.h"
#include "waves.h"
#include "utils.h"

int main(int argc, char * argv[]) {
	std::cout << "# Interstellar FErmi Diffusive Acceleration" << std::endl;

	if (argc == 2) {
//log_startup_information();

		Params params;
//params.set_from_file(argv[1]);
		params.set_params("u", 20.);
		params.set_params("v_A", 10.);
		params.set_params("dt", 1.);
		params.print();

		OutputManager O("output/test_skilling");

		Axis T(cgs::GeV);
		T.buildLogAxis(params.T_min, params.T_max, params.T_size);
		std::cout << " ... built energy axis T : " << T << "\n";

		Axis z(cgs::kpc);
		z.buildLinAxis(-params.H, params.H, params.z_size);
		std::cout << " ... built spatial axis z : " << z << "\n";

		Axis k(1. / cgs::pc);
		k.buildLogAxis(params.k_min, params.k_max, params.k_size);
		std::cout << " ... built wavenumber axis k : " << k << "\n";

		WaveSpectrum W(k, z, cgs::pc);
		W.build_initial_condition(params);
		std::cout << " ... built wave spectrum W : " << W << "\n";

		AdvectionVelocity u_z(z, cgs::km / cgs::sec);
		u_z.build_velocity_profile(params);
		std::cout << " ... built advection velocity u : " << u_z << "\n";

		DiffusionCoefficient D_zz(T, z, 1e28 * cgs::cm_2 / cgs::sec);
		std::cout << " ... built diffusion coefficient D_zz : " << D_zz << "\n";

		MomentumDiffusionCoefficient D_pp(T, z, cgs::GeV_c_2 / cgs::sec);
		std::cout << " ... built momentum diffusion coefficient D_pp : " << D_pp << "\n";

		EnergyLosses dpdt(T, z, cgs::GeV_c / cgs::sec);
		std::cout << " ... built energy losses dEdt : " << dpdt << "\n";

		Sources Q(T, z, 1. / cgs::cm_3 / cgs::GeV_c_3 / cgs::sec);
		std::cout << " ... built source term Q : " << Q << "\n";

		TimeEvolutor TEvolve(T, z, params.dt, 500 * 1000);
		TEvolve.doEvolveInSpace();
		TEvolve.doEvolveInMomentum();

		std::vector<Particle> particles;

		InjectionParams H_inj(0.1, 4.4);
		{
			PID pid = H1;
			std::cout << "running : " << pid << "\n";
			auto p = build_momentum_axis(T, pid);
			std::cout << "p : " << p << "\n";
			D_zz.build_Skilling(pid, params, W);
			std::cout << "D_zz : " << D_zz << "\n";
			D_pp.build_Skilling(pid, params, W);
			std::cout << "D_pp : " << D_pp << "\n";
			Q.build(pid, params, H_inj);
			std::cout << "Q : " << Q << "\n";
			dpdt.build(pid, params);
			std::cout << "dpdt : " << dpdt << "\n";
			O.dumpTimescales(Q, D_zz, D_pp, dpdt, u_z);
			O.dumpSpectrum(Q, D_zz, dpdt);
			//O.dumpProfile(Q, D_zz, dpdt, u_z);
			////Particle f_H = TEvolve.evolve(Q, D_zz, D_pp, dpdt, u_z, pid);
			////O.dumpSolution(f_H);
			//particles.emplace_back(f_H);
		}
	} else {
		std::cout << "Usage: ./idefa params.ini\n";
	}
	return 0;
}
