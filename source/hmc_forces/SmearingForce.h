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

namespace Update {

class SmearingForce : public StoutSmearing {
public:
	SmearingForce();
	~SmearingForce();

	void derivative(const extended_force_lattice_t& smearedDerivative, extended_force_lattice_t& unsmearedDerivative, extended_gauge_lattice_t& unsmearedLattice, real_t rho);
private:
	ForceVector ridder(extended_fermion_lattice_t& unsmearedLattice, int sited, unsigned int mud, int color, int site, unsigned int mu, real_t rho, real_t h = 0.05);

};

} /* namespace Update */
#endif /* SMEARINGFORCE_H_ */
