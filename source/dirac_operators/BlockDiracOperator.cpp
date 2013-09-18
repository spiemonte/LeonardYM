/*
 * BlockDiracOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "BlockDiracOperator.h"
#include "BlockDiracWilsonOperator.h"
#include "SquareBlockDiracWilsonOperator.h"

namespace Update {

BlockDiracOperator::BlockDiracOperator() : DiracOperator(), blockSize(4)  { }

BlockDiracOperator::BlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : DiracOperator(_lattice, _kappa), blockSize(4) { }

BlockDiracOperator::~BlockDiracOperator() { }

DiracOperator* BlockDiracOperator::getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters) {
	if (power == 1) {
		if (name == "DiracWilson") {
			BlockDiracWilsonOperator* result = new BlockDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			return result;
		}
		else {
			std::cout << "Dirac Wilson Operator" << name << " not supported or not supported in the block version!" << std::endl;
			exit(1);
		}
	} else if (power == 2) {
		if (name == "DiracWilson") {
			SquareBlockDiracWilsonOperator* result = new SquareBlockDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			return result;
		}
		else {
			std::cout << "Dirac Wilson Operator" << name << " not supported or not supported in the block version!" << std::endl;
			exit(1);
		}
	} else {
		std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
		exit(1);
	}
}

void BlockDiracOperator::setBlockSize(int _blockSize) {
	blockSize = _blockSize;
}

int BlockDiracOperator::getBlockSize() const {
	return blockSize;
}

} /* namespace Update */
