/*
 * FermionForce.h
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#ifndef FERMIONFORCE_H_
#define FERMIONFORCE_H_
#include "Environment.h"
#include "utils/LieGenerators.h"
#include "utils/ExpMap.h"

namespace Update {

class FermionForce {
public:
	FermionForce(real_t _kappa);
	virtual ~FermionForce();

	//Update only the fermion force, initialization to be called outside!
	virtual void derivative(extended_fermion_force_lattice_t& fermionForce, const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, real_t weight);

	virtual FermionicForceMatrix derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const = 0;

	GaugeGroup force(const environment_t& env, const FermionicForceMatrix& derivative, int site, unsigned int mu);

	virtual void setLattice(const extended_fermion_lattice_t& ) { }

protected:
	real_t kappa;
	
	ExponentialMap expMap;

	FermionicForceMatrix tensor(const GaugeVector& x, const GaugeVector& y) const;

	LieGenerator<FermionicGroup> fermionLieGenerator;
	LieGenerator<GaugeGroup> gaugeLieGenerator;
};

} /* namespace Update */
#endif /* FERMIONFORCE_H_ */
