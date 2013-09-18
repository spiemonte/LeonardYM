/*
 * SquareImprovedDiracWilsonOperator.cpp
 *
 *  Created on: Jul 16, 2012
 *      Author: spiem_01
 */

#include "SquareImprovedDiracWilsonOperator.h"
#include "ImprovedDiracWilsonOperator.h"
#include "../ImprovedFermionForce.h"

namespace Update {

SquareImprovedDiracWilsonOperator::SquareImprovedDiracWilsonOperator() : DiracOperator(), improvedDiracWilsonOperator(), csw(0) { }

SquareImprovedDiracWilsonOperator::SquareImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, real_t _csw, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), improvedDiracWilsonOperator(_lattice, _kappa, _csw, _gamma5), csw(_csw) { }

SquareImprovedDiracWilsonOperator::~SquareImprovedDiracWilsonOperator() { }

void SquareImprovedDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	improvedDiracWilsonOperator.setGamma5(gamma5);
	improvedDiracWilsonOperator.multiply(tmp, input);
	improvedDiracWilsonOperator.multiply(output, tmp);
}

void SquareImprovedDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	improvedDiracWilsonOperator.setGamma5(gamma5);
	improvedDiracWilsonOperator.multiply(tmp, vector1);
	improvedDiracWilsonOperator.multiplyAdd(output, tmp, vector2, alpha);
}

void SquareImprovedDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;//TODO
	improvedDiracWilsonOperator.setLattice(_lattice);
}

FermionForce* SquareImprovedDiracWilsonOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareImprovedDiracWilsonOperator, returning that for ImprovedDiracWilsonOperator" << std::endl;
	return new ImprovedFermionForce(kappa, csw);
}

void SquareImprovedDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	improvedDiracWilsonOperator.setKappa(_kappa);
}

real_t SquareImprovedDiracWilsonOperator::getCSW() const {
	return csw;
}

void SquareImprovedDiracWilsonOperator::setCSW(real_t _csw) {
	csw = _csw;
	improvedDiracWilsonOperator.setCSW(_csw);
}

} /* namespace Update */
