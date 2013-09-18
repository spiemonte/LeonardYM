/*
 * ImprovedFermionForce.h
 *
 *  Created on: Jun 11, 2012
 *      Author: spiem_01
 */

#ifndef IMPROVEDFERMIONFORCE_H_
#define IMPROVEDFERMIONFORCE_H_

#include "DiracWilsonFermionForce.h"

namespace Update {

class ImprovedFermionForce : public DiracWilsonFermionForce {
public:
	ImprovedFermionForce(real_t _kappa, real_t _csw = 1.);
	~ImprovedFermionForce();

	virtual FermionicForceMatrix derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const;
private:
	real_t csw;
};

} /* namespace Update */
#endif /* IMPROVEDFERMIONFORCE_H_ */
