/*
 * DiracWilsonOperator.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "TwistedDiracOperator.h"

namespace Update {

TwistedDiracOperator::TwistedDiracOperator() : DiracOperator(), diracOperator(0) { }

TwistedDiracOperator::TwistedDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

TwistedDiracOperator::~TwistedDiracOperator() { }

void TwistedDiracOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	diracOperator->multiply(output, input);
#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 2; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				output[site][mu][c] += std::complex<real_t>(0,twist)*input[site][mu][c];
			}
		}
		for (unsigned int mu = 2; mu < 4; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				output[site][mu][c] -= std::complex<real_t>(0,twist)*input[site][mu][c];
			}
		}
	}
}

void TwistedDiracOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	diracOperator->multiplyAdd(output, vector1, vector2, alpha);
#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 2; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				output[site][mu][c] += std::complex<real_t>(0,twist)*vector1[site][mu][c];
			}
		}
		for (unsigned int mu = 2; mu < 4; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				output[site][mu][c] -= std::complex<real_t>(0,twist)*vector1[site][mu][c];
			}
		}
	}
}

void TwistedDiracOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracOperator->setKappa(_kappa);
}

void TwistedDiracOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	diracOperator->setLattice(_lattice);
}

void TwistedDiracOperator::setDiracOperator(DiracOperator* op) {
	kappa = op->getKappa();
	//lattice = op->getLattice(); TODO
	diracOperator = op;
}

DiracOperator* TwistedDiracOperator::getDiracOperator() const {
	return diracOperator;
}

void TwistedDiracOperator::setTwist(real_t mu) {
	twist = mu;
}

real_t TwistedDiracOperator::getTwist() const {
	return twist;
}

FermionForce* TwistedDiracOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareTwistedDiracOperator, return that for DiracWilsonOperator" << std::endl;
	return diracOperator->getForce();
}

} /* namespace Update */
