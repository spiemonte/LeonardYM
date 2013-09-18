/*
 * Plaquette.cpp
 *
 *  Created on: Feb 29, 2012
 *      Author: spiem_01
 */

#include "Plaquette.h"
#include "StoutSmearing.h"
#include <omp.h>
#include <iostream>
#include "GlobalOutput.h"

namespace Update {

Plaquette::Plaquette() : LatticeSweep() { }

Plaquette::~Plaquette() { }

void Plaquette::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0.;

#pragma omp parallel for reduction(+:plaquette)
	for (unsigned int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaquette += real(trace(environment.gaugeLinkConfiguration[site][mu]*environment.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(environment.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(environment.gaugeLinkConfiguration[site][nu])));
			}
		}
	}
	reduceAllSum(plaquette);

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("plaquette");

		typedef extended_gauge_lattice_t::Layout Layout;
		std::cout << "Plaquette expectation value: " << plaquette/(6.*numberColors*Layout::globalVolume) << std::endl; //TODO
		output->write("plaquette", plaquette/(6.*numberColors*Layout::globalVolume));

		output->pop("plaquette");
	}

	/*StoutSmearing stout;
	fundamental_lattice_t result;
	stout.smearing(environment.gaugeLinkConfiguration, result, 0.15);
	std::cout << "Comparison smearing: " << this->temporalPlaquette(environment.gaugeLinkConfiguration) << " - " << this->temporalPlaquette(result) << std::endl;*/
}

long_real_t Plaquette::temporalPlaquette(const extended_gauge_lattice_t& gaugeLinkConfiguration) {
	typedef extended_gauge_lattice_t LT;
	long_real_t plaquette = 0.;
#pragma omp parallel for reduction(+:plaquette)
	for (unsigned int site = 0; site < gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 3; ++mu) {
			unsigned int nu = 3;
			plaquette += real(trace(gaugeLinkConfiguration[site][mu]*gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(gaugeLinkConfiguration[site][nu])));
		}
	}
	reduceAllSum(plaquette);
	if (isOutputProcess()) {
		typedef extended_gauge_lattice_t::Layout Layout;
		plaquette = plaquette/(3.*numberColors*Layout::globalVolume);
		std::cout << "Temporal Plaquette expectation value: " << plaquette << std::endl; //TODO
	}
	return plaquette; //TODO
}

} /* namespace Update */
