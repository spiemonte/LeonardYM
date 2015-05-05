/*
 * AdjointPolyakovLoop.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: spiem_01
 */

#include "AdjointPolyakovLoop.h"
#include "ConvertLattice.h"
#include "GlobalOutput.h"

namespace Update {

AdjointPolyakovLoop::AdjointPolyakovLoop() { }

AdjointPolyakovLoop::~AdjointPolyakovLoop() { }

void AdjointPolyakovLoop::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;
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

	for (int t = 0; t < Layout::glob_t; ++t) {
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndexT(site) == 0) {
				polyakov[site][3] = polyakov[site][3]*tmp[site][3];
			}
		}

		//Antialias
		swap = tmp;

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmp[site][mu] = swap[LT::sup(site,3)][mu];
			}
		}
		tmp.updateHalo();
	}

#pragma omp parallel for reduction(+:adjointPolyakovLoop)
	for (int site = 0; site < Layout::localsize; ++site) {
		if (Layout::globalIndexT(site) == 0) {
			AdjointGroup tmp;
			ConvertLattice<AdjointGroup,GaugeGroup>::toAdjoint(polyakov[site][3], tmp);
			adjointPolyakovLoop += trace(tmp);
		}
	}

	reduceAllSum(adjointPolyakovLoop);

	unsigned int spatialVolume = Layout::glob_spatial_volume;

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("polyakov_adjoint");

		std::cout << "Adjoint Polyakov Loop is (re) " << adjointPolyakovLoop/((numberColors*numberColors - 1)*spatialVolume) << std::endl;

		output->write("polyakov_adjoint", adjointPolyakovLoop/((numberColors*numberColors - 1)*spatialVolume));

		output->pop("polyakov_adjoint");
	}
}

} /* namespace Update */
