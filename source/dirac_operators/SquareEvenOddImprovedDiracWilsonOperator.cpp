/*
 * SquareImprovedDiracWilsonOperator.cpp
 *
 *  Created on: Jul 16, 2012
 *      Author: spiem_01
 */

#include "SquareEvenOddImprovedDiracWilsonOperator.h"
#include "ImprovedDiracWilsonOperator.h"
#include "hmc_forces/ImprovedFermionForce.h"

namespace Update {

SquareEvenOddImprovedDiracWilsonOperator::SquareEvenOddImprovedDiracWilsonOperator() : DiracOperator(), improvedDiracWilsonOperator(), csw(0) { }

SquareEvenOddImprovedDiracWilsonOperator::SquareEvenOddImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, real_t _csw, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), improvedDiracWilsonOperator(_lattice, _kappa, _csw, _gamma5), csw(_csw) { }

SquareEvenOddImprovedDiracWilsonOperator::~SquareEvenOddImprovedDiracWilsonOperator() { }

void SquareEvenOddImprovedDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	improvedDiracWilsonOperator.setGamma5(gamma5);
	improvedDiracWilsonOperator.multiply(tmp, input);
	improvedDiracWilsonOperator.multiply(output, tmp);
}

void SquareEvenOddImprovedDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	improvedDiracWilsonOperator.setGamma5(gamma5);
	improvedDiracWilsonOperator.multiply(tmp, vector1);
	improvedDiracWilsonOperator.multiply(output, tmp);
#pragma omp parallel for
	for (int site = 0; site < tmp.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] += alpha*vector2[site][mu];
	}
}

void SquareEvenOddImprovedDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;//TODO
	improvedDiracWilsonOperator.setLattice(_lattice);
}

FermionForce* SquareEvenOddImprovedDiracWilsonOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareEvenOddImprovedDiracWilsonOperator, returning that for ImprovedDiracWilsonOperator" << std::endl;
	return new ImprovedFermionForce(kappa, csw);
}

void SquareEvenOddImprovedDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	improvedDiracWilsonOperator.setKappa(_kappa);
}

real_t SquareEvenOddImprovedDiracWilsonOperator::getCSW() const {
	return csw;
}

void SquareEvenOddImprovedDiracWilsonOperator::setCSW(real_t _csw) {
	csw = _csw;
	improvedDiracWilsonOperator.setCSW(_csw);
}

} /* namespace Update */
