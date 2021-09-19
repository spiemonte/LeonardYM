#include "PolyakovLoop.h"
#include "io/GlobalOutput.h"

namespace Update {

PolyakovLoop::PolyakovLoop() { }

PolyakovLoop::~PolyakovLoop() { }

void PolyakovLoop::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	if (environment.measurement && isOutputProcess()) {
        GlobalOutput* output = GlobalOutput::getInstance();
        output->push("polyakov");
    }

	for (int nu = 3; nu >= 0; --nu) {
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

#pragma omp parallel for reduction(+:polyakovLoopRe,polyakovLoopIm)
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndex(site,nu) == 0) {
				std::complex<real_t> polyakovLoop = trace(polyakov[site][nu]);
				polyakovLoopRe += real(polyakovLoop);
				polyakovLoopIm += imag(polyakovLoop);
			}
		}

		reduceAllSum(polyakovLoopRe);
		reduceAllSum(polyakovLoopIm);

		unsigned int spatialVolume = Layout::glob_spatial_volume;

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			std::cout << "Polyakov Loop in the direction " << nu << " is (re) " << polyakovLoopRe/(numberColors*spatialVolume) << std::endl;
			std::cout << "Polyakov Loop in the direction " << nu << " is (im) " << polyakovLoopIm/(numberColors*spatialVolume) << std::endl;

			output->write("polyakov", polyakovLoopRe/(numberColors*spatialVolume));
			output->write("polyakov", polyakovLoopIm/(numberColors*spatialVolume));

		}
	}

	if (environment.measurement && isOutputProcess()) {
        GlobalOutput* output = GlobalOutput::getInstance();
        output->pop("polyakov");
	}
}

} /* namespace Update */
