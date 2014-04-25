/*
 * SquareDiracWilsonOperator.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "SquareTwistedDiracOperator.h"

namespace Update {

SquareTwistedDiracOperator::SquareTwistedDiracOperator() : DiracOperator(), diracOperator(0) { }

SquareTwistedDiracOperator::SquareTwistedDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

SquareTwistedDiracOperator::~SquareTwistedDiracOperator() { }

void SquareTwistedDiracOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	diracOperator->setGamma5(gamma5);
	diracOperator->multiply(output, input);
#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				output[site][mu][c] += twist*input[site][mu][c];
			}
		}
	}
}

void SquareTwistedDiracOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	diracOperator->setGamma5(gamma5);
	diracOperator->multiplyAdd(output, vector1, vector2, alpha + twist);
}

void SquareTwistedDiracOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracOperator->setKappa(_kappa);
}

void SquareTwistedDiracOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	diracOperator->setLattice(_lattice);
}

void SquareTwistedDiracOperator::setDiracOperator(DiracOperator* op) {
	kappa = op->getKappa();
	//lattice = op->getLattice(); TODO
	diracOperator = op;
}

DiracOperator* SquareTwistedDiracOperator::getDiracOperator() const {
	return diracOperator;
}

void SquareTwistedDiracOperator::setTwist(real_t mu) {
	twist = mu;
}

real_t SquareTwistedDiracOperator::getTwist() const {
	return twist;
}

FermionForce* SquareTwistedDiracOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareTwistedDiracOperator, return that for DiracWilsonOperator" << std::endl;
	return diracOperator->getForce();
}

} /* namespace Update */
