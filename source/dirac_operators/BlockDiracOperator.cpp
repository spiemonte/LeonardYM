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

BlockDiracOperator::BlockDiracOperator() : DiracOperator(), xBlockSize(4), yBlockSize(4), zBlockSize(4), tBlockSize(4)  { }

BlockDiracOperator::BlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : DiracOperator(_lattice, _kappa), xBlockSize(4), yBlockSize(4), zBlockSize(4), tBlockSize(4) { }

BlockDiracOperator::~BlockDiracOperator() { }

DiracOperator* BlockDiracOperator::getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters) {
	std::vector<unsigned int> blockSize = parameters.get< std::vector<unsigned int> >("deflation_block_size");
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

void BlockDiracOperator::setBlockSize(const std::vector<unsigned int>& _blockSize) {
	xBlockSize = _blockSize[0];
	yBlockSize = _blockSize[1];
	zBlockSize = _blockSize[2];
	tBlockSize = _blockSize[3];
}

std::vector<unsigned int> BlockDiracOperator::getBlockSize() const {
	std::vector<unsigned int> result;
	result.push_back(xBlockSize);
	result.push_back(yBlockSize);
	result.push_back(zBlockSize);
	result.push_back(tBlockSize);
	return result;
}

} /* namespace Update */
