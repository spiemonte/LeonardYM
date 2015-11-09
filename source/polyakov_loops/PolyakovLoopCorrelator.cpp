/*
 * PolyakovLoopCorrelator.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "PolyakovLoopCorrelator.h"
#include "io/GlobalOutput.h"
#include "utils/StoutSmearing.h"

namespace Update {

PolyakovLoopCorrelator::PolyakovLoopCorrelator() { }

PolyakovLoopCorrelator::~PolyakovLoopCorrelator() { }

void PolyakovLoopCorrelator::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;
	long_real_t polyakovLoopCorrelatorRe = 0;
	long_real_t polyakovLoopCorrelatorIm = 0;

	extended_gauge_lattice_t tmp = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t swap;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("level_stout_smearing_polyakov_correlator");
		double smearingRho = environment.configurations.get<double>("rho_stout_smearing");

		StoutSmearing stoutSmearing;
		for (unsigned int level = 0; level < numberLevelSmearing; ++level) {
			swap = tmp;
			stoutSmearing.smearing(swap, tmp, smearingRho);
		}
	} catch (NotFoundOption& ex) {
		if (isOutputProcess()) std::cout << "PolyakovLoopCorrelator::No smearing options found, proceeding without!" << std::endl;
	}
	
	extended_gauge_lattice_t polyakov;
	extended_gauge_lattice_t polyakov_translated;

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

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("polyakov_correlator");
	}
	
	polyakov_translated = polyakov;
	
	for (int l = 0; l < Layout::glob_x/2; ++l) {
		polyakovLoopCorrelatorRe = 0;
		polyakovLoopCorrelatorIm = 0;

		swap = polyakov_translated;
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				polyakov_translated[site][mu] = swap[LT::sup(site,0)][mu];
			}
		}
		polyakov_translated.updateHalo();
		
#pragma omp parallel for reduction(+:polyakovLoopCorrelatorRe,polyakovLoopCorrelatorIm)
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndexT(site) == 0) {
				std::complex<real_t> polyakovLoopCorrelator = conj(trace(polyakov[site][3]))*trace(polyakov_translated[site][3]);
				polyakovLoopCorrelatorRe += real(polyakovLoopCorrelator);
				polyakovLoopCorrelatorIm += imag(polyakovLoopCorrelator);
			}
		}

		reduceAllSum(polyakovLoopCorrelatorRe);
		reduceAllSum(polyakovLoopCorrelatorIm);

		unsigned int spatialVolume = Layout::glob_spatial_volume;

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			output->push("polyakov_correlator");

			std::cout << "PolyakovLoopCorrelator::Correlator at " << l << " " << polyakovLoopCorrelatorRe/(numberColors*spatialVolume) << " +I*" << polyakovLoopCorrelatorIm/(numberColors*spatialVolume) << std::endl;

			output->write("polyakov_correlator", polyakovLoopCorrelatorRe/(numberColors*spatialVolume));
			output->write("polyakov_correlator", polyakovLoopCorrelatorIm/(numberColors*spatialVolume));
			
			output->pop("polyakov_correlator");
		}
	
	}
	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->pop("polyakov_correlator");
	}
}

} /* namespace Update */
