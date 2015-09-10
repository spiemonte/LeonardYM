#include "Environment.h"
#include "utils/StoutSmearing.h"

namespace Update {

void Environment::synchronize() {
	try {
		int levels = configurations.get<int>("stout_smearing_levels");
		real_t rho = configurations.get<real_t>("stout_smearing_rho");
		if (levels > 1) {
			if (isOutputProcess()) std::cout << "Warning, number of levels " << levels << " not supported in unquenched mode!" << std::endl;
		}
		else if (levels == 0 || rho < 0.00000000001) {
			if (isOutputProcess()) std::cout << "Warning, comment smearing options if you don't use it!" << std::endl;
		}
		
	
		extended_gauge_lattice_t tmp = gaugeLinkConfiguration, smeared;
		StoutSmearing stoutSmearing;
	
		for (int level = 0; level < levels; ++level) {
			stoutSmearing.smearing(tmp, smeared, rho);
			tmp = smeared;
		}
	
#ifdef ADJOINT
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::convert(fermionicLinkConfiguration, smeared);
#endif
#ifndef ADJOINT
		fermionicLinkConfiguration = smeared;
#endif
	}
	catch (NotFoundOption& ex) {
#ifdef ADJOINT
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::convert(fermionicLinkConfiguration, gaugeLinkConfiguration);
#endif
#ifndef ADJOINT
		fermionicLinkConfiguration = gaugeLinkConfiguration;
#endif
	}

	try {
		std::string bc = configurations.get<std::string>("boundary_conditions");
		if (bc == "fermion_antiperiodic") {
			switchAntiperiodicBc(fermionicLinkConfiguration);
		}
		else if (bc == "fermion_spatial_antiperiodic") {
			switchSpatialAntiperiodicBc(fermionicLinkConfiguration);
		}
		else if (bc == "fermion_periodic") {
		}
		else if (bc == "open") {
			switchOpenBc(fermionicLinkConfiguration);
			switchOpenBc(gaugeLinkConfiguration);
		}
		else {
			static int count = 0;
			if (count == 0 && isOutputProcess()) {
				std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
				count = count + 1;
			}
			switchAntiperiodicBc(fermionicLinkConfiguration);
		}
	}
	catch (NotFoundOption& ex) {
		static int count = 0;
		if (count == 0 && isOutputProcess()) {
			std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
			count = count + 1;
		}
		switchAntiperiodicBc(fermionicLinkConfiguration);
	}
}

}
