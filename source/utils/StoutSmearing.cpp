/*
 * StoutSmearing.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#include "StoutSmearing.h"
#include "actions/WilsonGaugeAction.h"
#include "utils/ToString.h"

namespace Update {

StoutSmearing::StoutSmearing() : ExponentialMap(), wga(0.) { }

StoutSmearing::~StoutSmearing() { }

void StoutSmearing::smearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, real_t rho) {
#pragma omp parallel for
	for (int site = 0; site < input.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			output[site][mu] = this->smearLink(input, site, mu, rho);
		}
	}
	output.updateHalo();
}

void StoutSmearing::spatialSmearing(const extended_gauge_lattice_t& initialInput, extended_gauge_lattice_t& output, unsigned int numberLevels, real_t rho) {
	extended_gauge_lattice_t input = initialInput;
	typedef extended_gauge_lattice_t LT;
	for (unsigned int level = 0; level < numberLevels; ++level) {
#pragma omp parallel for
		for (int site = 0; site < input.localsize; ++site) {
			for (unsigned int mu = 0; mu < 3; ++mu) {
				GaugeGroup staple;
				set_to_zero(staple);
				for (unsigned int nu = 0; nu < 3; ++nu) {
					if (nu != mu) {
						staple += input[LT::sup(site,mu)][nu]*htrans(input[LT::sup(site,nu)][mu])*htrans(input[site][nu]);
						staple += htrans(input[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(input[LT::sdn(site,nu)][mu])*input[LT::sdn(site,nu)][nu];
					}
				}
				GaugeGroup omega = rho*staple*htrans(input[site][mu]);
				GaugeGroup iStout = 0.5*(htrans(omega) - omega);
				GaugeGroup toExp = -(iStout - (trace(iStout)/static_cast<real_t>(numberColors))*identity);
				output[site][mu] = this->exp(toExp)*input[site][mu];
			}
			output[site][3] = input[site][3];
		}
		output.updateHalo();
		input = output;
	}
}

GaugeGroup StoutSmearing::smearLink(const extended_gauge_lattice_t& input, int site, unsigned int mu, real_t rho) const {
	GaugeGroup staple = htrans(wga.staple(input,site,mu));
	GaugeGroup omega = rho*staple*htrans(input[site][mu]);
	omega = rho*input[site][mu];//TODO
	GaugeGroup iStout = 0.5*(htrans(omega) - omega);

	GaugeGroup toExp = -(iStout - (trace(iStout)/static_cast<real_t>(numberColors))*identity);

	return input[site][mu];//*input[site][mu];//this->exp(toExp)*input[site][mu];
}

#ifdef ADJOINT
FermionicGroup StoutSmearing::smearLink(const extended_fermion_lattice_t& input, int site, unsigned int mu, real_t rho) const {
	typedef extended_fermion_lattice_t LT;
	FermionicGroup staple;
	set_to_zero(staple);
	for (unsigned int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			staple += input[LT::sup(site,mu)][nu]*htrans(input[LT::sup(site,nu)][mu])*htrans(input[site][nu]);
			staple += htrans(input[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(input[LT::sdn(site,nu)][mu])*input[LT::sdn(site,nu)][nu];
		}
	}
	FermionicGroup omega = rho*staple*htrans(input[site][mu]);
	FermionicGroup iStout = 0.5*(htrans(omega) - omega);
	FermionicGroup toExp = -(iStout - (trace(iStout)/static_cast<real_t>(numberColors))*adjoint_identity);
	return /*input[site][mu]*/input[site][mu];//(this->exp(toExp)*input[site][mu]);
}
#endif

} /* namespace Update */

