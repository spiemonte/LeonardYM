/*
 * BasicSquareDiracWilsonOperator.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: spiem_01
 */

#include "BasicSquareDiracWilsonOperator.h"
#include "BasicDiracWilsonOperator.h"
#include "../DiracWilsonFermionForce.h"

namespace Update {

BasicSquareDiracWilsonOperator::BasicSquareDiracWilsonOperator() : DiracOperator() { }

BasicSquareDiracWilsonOperator::BasicSquareDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

BasicSquareDiracWilsonOperator::~BasicSquareDiracWilsonOperator() { }

void BasicSquareDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	BasicDiracWilsonOperator diracWilsonOperator(lattice, kappa, gamma5);
	diracWilsonOperator.multiply(tmp, input);
	diracWilsonOperator.multiply(output, tmp);
}

void BasicSquareDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	BasicDiracWilsonOperator diracWilsonOperator(lattice, kappa, gamma5);
	diracWilsonOperator.multiply(tmp, vector1);
	diracWilsonOperator.multiplyAdd(output, tmp, vector2, alpha);
}

FermionForce* BasicSquareDiracWilsonOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareDiracWilsonOperator, return that for DiracWilsonOperator" << std::endl;
	return new DiracWilsonFermionForce(kappa);
}

}
