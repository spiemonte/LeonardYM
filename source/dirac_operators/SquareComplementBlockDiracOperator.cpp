/*
 * SquareComplementBlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "SquareComplementBlockDiracOperator.h"
#include "../DiracWilsonFermionForce.h"
#include "../AlgebraUtils.h"

namespace Update {

SquareComplementBlockDiracOperator::SquareComplementBlockDiracOperator(ComplementBlockDiracOperator* _K) : BlockDiracOperator(), K(_K) { }

SquareComplementBlockDiracOperator::~SquareComplementBlockDiracOperator() { }

void SquareComplementBlockDiracOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
#pragma omp parallel for
	for (int site = 0; site < input.completesize; ++site) {
		tmpVector[site][0] = input[site][0];
		tmpVector[site][1] = input[site][1];
		tmpVector[site][2] = -input[site][2];
		tmpVector[site][3] = -input[site][3];
	}

	K->multiply(output,tmpVector);

#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		tmpVector[site][0] = output[site][0];
		tmpVector[site][1] = output[site][1];
		tmpVector[site][2] = -output[site][2];
		tmpVector[site][3] = -output[site][3];
	}

	K->multiply(output,tmpVector);
}

void SquareComplementBlockDiracOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
#pragma omp parallel for
	for (int site = 0; site < vector1.completesize; ++site) {
		tmpVector[site][0] = vector1[site][0];
		tmpVector[site][1] = vector1[site][1];
		tmpVector[site][2] = -vector1[site][2];
		tmpVector[site][3] = -vector1[site][3];
	}

	K->multiply(output,tmpVector);

#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		tmpVector[site][0] = output[site][0];
		tmpVector[site][1] = output[site][1];
		tmpVector[site][2] = -output[site][2];
		tmpVector[site][3] = -output[site][3];
	}

	K->multiply(output,tmpVector);

	//Now we add the shift
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] += alpha*vector2[site][mu];
	}
}

void SquareComplementBlockDiracOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	K->setLattice(_lattice);
}

FermionForce* SquareComplementBlockDiracOperator::getForce() const {
	return NULL;
}

void SquareComplementBlockDiracOperator::setPrecision(const real_t& _precision) {
	K->setPrecision(_precision);
}

void SquareComplementBlockDiracOperator::setMaximumSteps(int steps) {
	K->setMaximumSteps(steps);
}

void SquareComplementBlockDiracOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	K->setKappa(_kappa);
}

} /* namespace Update */
