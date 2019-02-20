/*
 * OverlapFermionForce.h
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#ifndef OVERLAPFERMIONFORCE_H_
#define OVERLAPFERMIONFORCE_H_

#include "dirac_operators/DiracWilsonOperator.h"
#include "dirac_operators/SquareDiracWilsonOperator.h"
#include "DiracWilsonFermionForce.h"
#include "dirac_functions/Polynomial.h"

namespace Update {

class OverlapFermionForce : public Update::DiracWilsonFermionForce {
public:
	OverlapFermionForce(real_t _kappa, real_t _mass, Polynomial const* _squareRootApproximation);
	~OverlapFermionForce();

	virtual void derivative(extended_fermion_force_lattice_t& fermionForce, const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, real_t weight);
private:
	real_t mass;
	DiracWilsonOperator* diracWilsonOperator;
	SquareDiracWilsonOperator* squareDiracWilsonOperator;
	Polynomial const* squareRootApproximation;

	std::vector<extended_dirac_vector_t> left_dirac_vectors;
	std::vector<extended_dirac_vector_t> right_dirac_vectors;
	extended_dirac_vector_t tmp;
};

} /* namespace Update */
#endif /* DIRACWILSONFERMIONFORCE_H_ */
