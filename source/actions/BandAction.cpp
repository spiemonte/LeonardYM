/*
 * BandAction.cpp
 *
 *  Created on: Jan 11, 2013
 *      Author: spiem_01
 */

#include "BandAction.h"

namespace Update {

BandAction::BandAction(Force* _subForce) : subForce(_subForce) { }

BandAction::~BandAction() { }

GaugeGroup BandAction::force(const environment_t& env, int site, int mu) const {
	return subForce->force(env, site, mu);
}

void BandAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	subForce->updateForce(forceLattice,env);
	typedef extended_gauge_lattice_t::Layout Layout;
	std::vector< std::pair<int, int> >::const_iterator i = bands.begin();
	
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		bool isActive = false;
		int tGlobalCoordinate = Layout::globalIndexT(site);
		for (i = bands.begin(); i != bands.end(); ++i) {
			if (tGlobalCoordinate >= i->first && tGlobalCoordinate < i->second) {
				isActive = true;
				break;
			}
		}
		if (!isActive) {
			for (unsigned int mu = 0; mu < 4; ++mu) set_to_zero(forceLattice[site][mu]);
		}
	}
}

void BandAction::addBand(int t0, int t1) {
	bands.push_back(std::pair<int, int>(t0,t1));
}


} /* namespace Update */
