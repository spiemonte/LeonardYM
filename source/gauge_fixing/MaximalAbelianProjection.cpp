#include "MaximalAbelianProjection.h"
#include "utils/LieGenerators.h"

namespace Update {

MaximalAbelianProjection::MaximalAbelianProjection() : LatticeSweep() { }

MaximalAbelianProjection::MaximalAbelianProjection(const MaximalAbelianProjection& toCopy) : LatticeSweep(toCopy) { }

MaximalAbelianProjection::~MaximalAbelianProjection() { }

void MaximalAbelianProjection::execute(environment_t& environment) {
#pragma omp parallel for
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int c = 0; c < numberColors; ++c) {
				for (unsigned int d = c+1; d < numberColors; ++d) {
					environment.gaugeLinkConfiguration[site][mu].at(c,d) = 0.;
					environment.gaugeLinkConfiguration[site][mu].at(d,c) = 0.;
				}
			}
		}
	}

	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
}

}
