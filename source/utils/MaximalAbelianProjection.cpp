#include "MaximalAbelianProjection.h"
#include "LieGenerators.h"

namespace Update {

MaximalAbelianProjection::MaximalAbelianProjection() : LatticeSweep() { }

MaximalAbelianProjection::MaximalAbelianProjection(const MaximalAbelianProjection& toCopy) : LatticeSweep(toCopy) { }

MaximalAbelianProjection::~MaximalAbelianProjection() { }

void MaximalAbelianProjection::execute(environment_t& environment) {
#if NUMCOLORS == 2
	
#pragma omp parallel for
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			environment.gaugeLinkConfiguration[site][mu].at(1,0) = 0.;
			environment.gaugeLinkConfiguration[site][mu].at(0,1) = 0.;
		}
	}

	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
#endif
}

void MaximalAbelianProjection::registerParameters(po::options_description&) { }

}

