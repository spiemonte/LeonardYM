/*
 * StoutSmearing.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#ifndef STOUTSMEARING_H_
#define STOUTSMEARING_H_
#include "Environment.h"
#include "actions/WilsonGaugeAction.h"
#include "ExpMap.h"

namespace Update {

class StoutSmearing : public ExponentialMap {
public:
	StoutSmearing();
	~StoutSmearing();

	void smearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, real_t rho);
	void spatialSmearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, unsigned int numberLevels, real_t rho);

protected:
	GaugeGroup smearLink(const extended_gauge_lattice_t& input, int site, unsigned int mu, real_t rho) const;
#ifdef ADJOINT
	FermionicGroup smearLink(const extended_fermion_lattice_t& input, int site, unsigned int mu, real_t rho) const;
#endif

private:
	WilsonGaugeAction wga;
};

} /* namespace Update */
#endif /* STOUTSMEARING_H_ */
