/*
 * SmearingForce.h
 *
 *  Created on: Dec 7, 2012
 *      Author: spiem_01
 */

#ifndef SMEARINGFORCE_H_
#define SMEARINGFORCE_H_
#include "utils/StoutSmearing.h"
#include "utils/EigenTraits.h"
#include "utils/ExpMap.h"

namespace Update {

class SmearingForce : public StoutSmearing {
public:
	SmearingForce();
	~SmearingForce();

	void force(const extended_fermion_force_lattice_t& actionDerivative, const extended_gauge_lattice_t& unsmearedLattice, extended_gauge_lattice_t& unsmearedDerivative, real_t rho);
private:
	real_t ridder(const extended_fermion_force_lattice_t& actionDerivative, extended_gauge_lattice_t& unsmearedLattice, int sited, unsigned int mud, int color, int site, unsigned int mu, real_t rho, real_t h = 0.2);

	ExponentialMap expMap;
	LieGenerator<GaugeGroup> gaugeLieGenerators;
};

} /* namespace Update */
#endif /* SMEARINGFORCE_H_ */
