/*
 * Force.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#include "Force.h"

namespace Update {

Force::Force() { }

Force::~Force() { }

void Force::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	//Calculate the force
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			forceLattice[site][mu] = this->force(env, site, mu);
		}
	}
	forceLattice.updateHalo();//TODO maybe not needed?
}

} /* namespace Update */
