/*
 * WilsonGaugeAction.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#include "WilsonGaugeAction.h"

namespace Update {

WilsonGaugeAction::WilsonGaugeAction(real_t _beta) : GaugeAction(_beta) { }

WilsonGaugeAction::~WilsonGaugeAction() { }

GaugeGroup WilsonGaugeAction::staple(const extended_gauge_lattice_t& lattice, int site, int mu) const {
	typedef extended_gauge_lattice_t LT;
	GaugeGroup ris;
	set_to_zero(ris);
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			ris += lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			ris += htrans(lattice[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(lattice[LT::sdn(site,nu)][mu])*lattice[LT::sdn(site,nu)][nu];
		}
	}
	return ris;
}

GaugeGroup WilsonGaugeAction::force(const extended_gauge_lattice_t& lattice, int site, int mu) const {
	GaugeGroup plaquette = lattice[site][mu]*this->staple(lattice, site, mu);
	GaugeGroup force = -(0.25*this->getBeta()/numberColors)*(htrans(plaquette) - plaquette);
	std::complex<real_t> trc = trace(force);
	//Traceless part
	for (int i = 0; i < numberColors; ++i) {
		force.at(i,i) -= std::complex<real_t>(real(trc)/numberColors,imag(trc)/numberColors);
	}
	return force;
}

long_real_t WilsonGaugeAction::energy(const environment_t& env) {
	typedef extended_gauge_lattice_t LT;
	long_real_t energy = 0.;
#pragma omp parallel for reduction(+:energy)
	for (int site = 0; site < env.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += env.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(env.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(env.gaugeLinkConfiguration[site][nu]);
			}
			energy += -(this->getBeta()/numberColors)*real(trace(env.gaugeLinkConfiguration[site][mu]*plaqs));
		}
	}
	reduceAllSum(energy);
	return energy;
}

real_t WilsonGaugeAction::deltaAction(const extended_gauge_lattice_t& lattice, const GaugeGroup& trial, const GaugeGroup& staple, int site, int mu) const {
	real_t oldAction = (real(trace(lattice[site][mu]*staple)));
	real_t newAction = (real(trace(trial*staple)));
	return -this->getBeta()*(newAction-oldAction)/(numberColors);
}

} /* namespace Update */
