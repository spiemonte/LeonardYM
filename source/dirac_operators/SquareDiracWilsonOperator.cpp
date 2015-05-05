/*
 * SquareDiracWilsonOperator.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "SquareDiracWilsonOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"

namespace Update {

SquareDiracWilsonOperator::SquareDiracWilsonOperator() : DiracOperator(), diracWilsonOperator() { }

SquareDiracWilsonOperator::SquareDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), diracWilsonOperator(_lattice, _kappa, _gamma5) { }

SquareDiracWilsonOperator::~SquareDiracWilsonOperator() { }

void SquareDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	diracWilsonOperator.setGamma5(gamma5);
	diracWilsonOperator.multiply(tmp, input);
	diracWilsonOperator.multiply(output, tmp);
}

void SquareDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	diracWilsonOperator.setGamma5(gamma5);
	diracWilsonOperator.multiply(tmp, vector1);
	diracWilsonOperator.multiplyAdd(output, tmp, vector2, alpha);
}

void SquareDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracWilsonOperator.setKappa(_kappa);
}

void SquareDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	diracWilsonOperator.setLattice(_lattice);
}

FermionForce* SquareDiracWilsonOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareDiracWilsonOperator, return that for DiracWilsonOperator" << std::endl;
	return new DiracWilsonFermionForce(kappa);
}

} /* namespace Update */
