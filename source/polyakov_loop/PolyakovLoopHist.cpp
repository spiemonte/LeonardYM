/*
 * PolyakovLoop.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "PolyakovLoop.h"
#include "GlobalOutput.h"

namespace Update {

PolyakovLoop::PolyakovLoop():hist2d_(NULL),evhist_(NULL),evhist2d_(NULL) { }

PolyakovLoop::~PolyakovLoop() {
    //add here histogrammer: collect from all proc and write out
    // delete of histogrammer
}

void PolyakovLoop::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;
	long_real_t polyakovLoopRe = 0;
	long_real_t polyakovLoopIm = 0;

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

	bool histogram(false);
	bool histogramev(false);
	try {
	        histogram = environment.configurations.get<bool >("histpgram_polyakov_loop");
            histogramev = environment.configurations.get<bool >("histpgram_polyakov_loop");
            if (isOutputProcess()) std::cout <<"Histograms of PlLoop prepared"<<std::endl;
	}catch(NotFoundOption& ex){
	    histogram =histogramev=false;
	}
	if(histogram){
	    if(!hist2d_){

	    }
	 // add here histogrammer
	}
	if(histogramev){
	 // add here histogrammer
	}

#pragma omp parallel for reduction(+:polyakovLoopRe,polyakovLoopIm)
	for (int site = 0; site < Layout::localsize; ++site) {
		if (Layout::globalIndexT(site) == 0) {
			std::complex<real_t> polyakovLoop = trace(polyakov[site][3]);
			polyakovLoopRe += real(polyakovLoop);
			polyakovLoopIm += imag(polyakovLoop);
		}
	}

	reduceAllSum(polyakovLoopRe);
	reduceAllSum(polyakovLoopIm);

	unsigned int spatialVolume = Layout::glob_spatial_volume;

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("polyakov");

		std::cout << "Polyakov Loop is (re) " << polyakovLoopRe/(numberColors*spatialVolume) << std::endl;
		std::cout << "Polyakov Loop is (im) " << polyakovLoopIm/(numberColors*spatialVolume) << std::endl;

		output->write("polyakov", polyakovLoopRe/(numberColors*spatialVolume));
		output->write("polyakov", polyakovLoopIm/(numberColors*spatialVolume));

		output->pop("polyakov");
	}
}

} /* namespace Update */
