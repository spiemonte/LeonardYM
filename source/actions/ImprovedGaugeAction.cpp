#include "ImprovedGaugeAction.h"

namespace Update {

ImprovedGaugeAction::ImprovedGaugeAction(real_t _beta, real_t _u0) : GaugeAction(_beta), u0(_u0)  { }

ImprovedGaugeAction::~ImprovedGaugeAction() { }

GaugeGroup ImprovedGaugeAction::staple(const extended_gauge_lattice_t& lattice, int site, int mu) const {
	typedef extended_gauge_lattice_t LT;
	GaugeGroup ris;
	set_to_zero(ris);
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			ris += (5./3.)*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			ris += (5./3.)*htrans(lattice[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(lattice[LT::sdn(site,nu)][mu])*lattice[LT::sdn(site,nu)][nu];

			ris -= (1./(12.*u0*u0))* (lattice[LT::sup(site, mu)][mu])*(lattice[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(lattice[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[site][nu]);
			ris -= (1./(12.*u0*u0))* (lattice[LT::sup(site, mu)][nu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[LT::sup(LT::sdn(site, mu), nu)][mu])*htrans(lattice[LT::sdn(site, mu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
			ris -= (1./(12.*u0*u0))* (lattice[LT::sup(site, mu)][mu])*htrans(lattice[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu])*htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][mu])*htrans(lattice[LT::sdn(site, nu)][mu])*(lattice[LT::sdn(site, nu)][nu]);
			ris -= (1./(12.*u0*u0))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sdn(site, nu)][mu])*htrans(lattice[LT::sdn(LT::sdn(site, mu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, mu), nu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
			ris -= (1./(12.*u0*u0))* (lattice[LT::sup(site, mu)][nu])*(lattice[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(lattice[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][nu])*htrans(lattice[site][nu]);
			ris -= (1./(12.*u0*u0))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu])*htrans(lattice[LT::sdn(LT::sdn(site, nu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, nu), nu)][nu])*(lattice[LT::sdn(site, nu)][nu]);
		}
	}
	return ris;
}

GaugeGroup ImprovedGaugeAction::force(const extended_gauge_lattice_t& lattice, int site, int mu) const {
	GaugeGroup plaquette = lattice[site][mu]*this->staple(lattice, site, mu);
	GaugeGroup force = -(0.25*this->getBeta()/numberColors)*(htrans(plaquette) - plaquette);
	std::complex<real_t> trc = trace(force);
	//Traceless part
	for (int i = 0; i < numberColors; ++i) {
		force.at(i,i) -= std::complex<real_t>(real(trc)/numberColors,imag(trc)/numberColors);
	}
	return force;
}

long_real_t ImprovedGaugeAction::energy(const environment_t& env) {
	typedef extended_gauge_lattice_t LT;
	long double energy = 0.;
#pragma omp parallel for reduction(+: energy)
	for (int site = 0; site < env.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += (5./3.)*env.gaugeLinkConfiguration[LT::sup(site,mu)][nu]*htrans(env.gaugeLinkConfiguration[LT::sup(site,nu)][mu])*htrans(env.gaugeLinkConfiguration[site][nu]);
				plaqs -= (1./(12.*u0*u0))*(env.gaugeLinkConfiguration[LT::sup(site, mu)][mu])*(env.gaugeLinkConfiguration[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(env.gaugeLinkConfiguration[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(env.gaugeLinkConfiguration[LT::sup(site, nu)][mu])*htrans(env.gaugeLinkConfiguration[site][nu]);
				plaqs -= (1./(12.*u0*u0))*(env.gaugeLinkConfiguration[LT::sup(site, mu)][nu])*(env.gaugeLinkConfiguration[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(env.gaugeLinkConfiguration[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(env.gaugeLinkConfiguration[LT::sup(site, nu)][nu])*htrans(env.gaugeLinkConfiguration[site][nu]);
			}
			energy += -(this->getBeta()/numberColors)*real(trace(env.gaugeLinkConfiguration[site][mu]*plaqs));
		}
	}
	reduceAllSum(energy);
	return energy;
}

real_t ImprovedGaugeAction::deltaAction(const extended_gauge_lattice_t& lattice, const GaugeGroup& trial, const GaugeGroup& staple, int site, int mu) const {
	real_t oldAction = (real(trace(lattice[site][mu]*staple)));
	real_t newAction = (real(trace(trial*staple)));
	return -this->getBeta()*(newAction-oldAction)/(numberColors);
}

} /* namespace Update */
