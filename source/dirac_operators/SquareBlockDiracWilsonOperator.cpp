/*
 * SquareBlockDiracWilsonOperator.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "SquareBlockDiracWilsonOperator.h"
#include "../BlockDiracWilsonFermionForce.h"

namespace Update {

SquareBlockDiracWilsonOperator::SquareBlockDiracWilsonOperator() : BlockDiracOperator(), blockDiracWilsonOperator() { }

SquareBlockDiracWilsonOperator::SquareBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : BlockDiracOperator(_lattice, _kappa), blockDiracWilsonOperator(_lattice, _kappa) { }

SquareBlockDiracWilsonOperator::~SquareBlockDiracWilsonOperator() { }

void SquareBlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	blockDiracWilsonOperator.multiply(tmp, input);
	blockDiracWilsonOperator.multiply(output, tmp);
}

void SquareBlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	blockDiracWilsonOperator.multiply(tmp, vector1);
	blockDiracWilsonOperator.multiplyAdd(output, tmp, vector2, alpha);
}

void SquareBlockDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	blockDiracWilsonOperator.setKappa(_kappa);
}

void SquareBlockDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	blockDiracWilsonOperator.setLattice(_lattice);
}

FermionForce* SquareBlockDiracWilsonOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareBlockDiracWilsonOperator, return that for DiracWilsonOperator" << std::endl;
	return new BlockDiracWilsonFermionForce(kappa, blockSize);
}

void SquareBlockDiracWilsonOperator::setBlockSize(int _blockSize) {
	blockSize = _blockSize;
	blockDiracWilsonOperator.setBlockSize(_blockSize);
}

} /* namespace Update */
