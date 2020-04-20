#include "AdjointPolyakovLoop.h"
#include "utils/ConvertLattice.h"
#include "io/GlobalOutput.h"

namespace Update {

AdjointPolyakovLoop::AdjointPolyakovLoop() { }

AdjointPolyakovLoop::~AdjointPolyakovLoop() { }

void AdjointPolyakovLoop::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	if (environment.measurement && isOutputProcess()) {
                GlobalOutput* output = GlobalOutput::getInstance();
                output->push("polyakov_adjoint");
        }


	for (int nu = 3; nu >= 0; --nu) {
		long_real_t adjointPolyakovLoop = 0;

		extended_gauge_lattice_t tmp = environment.gaugeLinkConfiguration;
		extended_gauge_lattice_t swap;
		extended_gauge_lattice_t polyakov;

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_identity(polyakov[site][mu]);
			}
		}

		for (int t = 0; t < Layout::glob[nu]; ++t) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				if (Layout::globalIndex(site,nu) == 0) {
					polyakov[site][nu] = polyakov[site][nu]*tmp[site][nu];
				}
			}

			//Antialias
			swap = tmp;

#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					tmp[site][mu] = swap[LT::sup(site,nu)][mu];
				}
			}
			tmp.updateHalo();
		}

#pragma omp parallel for reduction(+:adjointPolyakovLoop)
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndex(site,nu) == 0) {
				AdjointGroup tmp;
				ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(polyakov[site][nu], tmp);
				adjointPolyakovLoop += trace(tmp);
			}
		}

		reduceAllSum(adjointPolyakovLoop);

		unsigned int spatialVolume = Layout::glob_spatial_volume;

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			std::cout << "Adjoint Polyakov Loop in the direction " << nu << " is (re) " << adjointPolyakovLoop/((numberColors*numberColors - 1)*spatialVolume) << std::endl;

			output->write("polyakov_adjoint", adjointPolyakovLoop/((numberColors*numberColors - 1)*spatialVolume));
		}
	}

	if (environment.measurement && isOutputProcess()) {
                GlobalOutput* output = GlobalOutput::getInstance();
                output->pop("polyakov_adjoint");
	}
}

} /* namespace Update */
