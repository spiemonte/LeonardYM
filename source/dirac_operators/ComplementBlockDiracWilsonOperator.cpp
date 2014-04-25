/*
 * ComplementBlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "ComplementBlockDiracWilsonOperator.h"
#include "../DiracWilsonFermionForce.h"
#include "../AlgebraUtils.h"

namespace Update {

ComplementBlockDiracWilsonOperator::ComplementBlockDiracWilsonOperator() : DiracOperator(), diracWilsonOperator(), blockDiracWilsonOperator(), squareBlockDiracWilsonOperator(), log(false), counterSteps(0) { }

ComplementBlockDiracWilsonOperator::ComplementBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : DiracOperator(_lattice, _kappa), diracWilsonOperator(_lattice, _kappa), blockDiracWilsonOperator(_lattice, _kappa), squareBlockDiracWilsonOperator(_lattice, _kappa), log(false), counterSteps(0) { }

ComplementBlockDiracWilsonOperator::~ComplementBlockDiracWilsonOperator() { }

void ComplementBlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	//First we apply D
	diracWilsonOperator.multiply(output, input);
	
	//Now we apply D_1^(-1) to (D.input)
	biConjugateGradient.solve(&squareBlockDiracWilsonOperator, output, tmpVector);
	blockDiracWilsonOperator.multiply(output, tmpVector);
	
	if (log && isOutputProcess()) std::cout << "ComplementBlockDiracWilsonOperator::Complement inner inversion done in: " << biConjugateGradient.getLastSteps() << std::endl;
}

void ComplementBlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	//First we apply D
	diracWilsonOperator.multiply(output, vector1);
	
	//Now we apply D_1^(-1) to (D.input)
	biConjugateGradient.solve(&squareBlockDiracWilsonOperator, output, tmpVector);
	counterSteps += biConjugateGradient.getLastSteps();	
	blockDiracWilsonOperator.multiply(output, tmpVector);

	//Now we add the shift
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] += alpha*vector2[site][mu];
	}
	
	if (log && isOutputProcess()) std::cout << "ComplementBlockDiracWilsonOperator::Complement inner inversion done in: " << biConjugateGradient.getLastSteps() << std::endl;
	
}

void ComplementBlockDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	blockDiracWilsonOperator.setLattice(_lattice);
	squareBlockDiracWilsonOperator.setLattice(_lattice);
	diracWilsonOperator.setLattice(_lattice);
}

void ComplementBlockDiracWilsonOperator::setBlockSize(const std::vector<unsigned int>& _blockSize) {
	blockDiracWilsonOperator.setBlockSize(_blockSize);
	squareBlockDiracWilsonOperator.setBlockSize(_blockSize);
}

std::vector<unsigned int> ComplementBlockDiracWilsonOperator::getBlockSize() const {
	return blockDiracWilsonOperator.getBlockSize();
}

FermionForce* ComplementBlockDiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

void ComplementBlockDiracWilsonOperator::setPrecision(const real_t& _precision) {
	biConjugateGradient.setPrecision(_precision);
}

real_t ComplementBlockDiracWilsonOperator::getPrecision() const {
	return biConjugateGradient.getPrecision();
}

void ComplementBlockDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	blockDiracWilsonOperator.setKappa(_kappa);
	diracWilsonOperator.setKappa(_kappa);
	squareBlockDiracWilsonOperator.setKappa(_kappa);
}

void ComplementBlockDiracWilsonOperator::setLog(bool _log) {
	log = _log;
}

void ComplementBlockDiracWilsonOperator::resetCounterInnerSteps() {
	counterSteps = 0;
}

int ComplementBlockDiracWilsonOperator::getInnerSteps() const {
	return counterSteps;
}

int ComplementBlockDiracWilsonOperator::getInnerInverterSteps() const {
	return biConjugateGradient.getLastSteps();
}


} /* namespace Update */
