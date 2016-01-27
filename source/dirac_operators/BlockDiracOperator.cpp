/*
 * BlockDiracOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "BlockDiracOperator.h"
#include "BlockDiracWilsonOperator.h"
#include "BlockImprovedDiracWilsonOperator.h"
#include "SquareBlockDiracWilsonOperator.h"

namespace Update {

BlockDiracOperator::BlockDiracOperator(Color _color) : DiracOperator(), xBlockSize(6), yBlockSize(6), zBlockSize(6), tBlockSize(6), color(_color) {
	this->initialize_index_lattice();
}

BlockDiracOperator::BlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, Color _color) : DiracOperator(_lattice, _kappa), xBlockSize(6), yBlockSize(6), zBlockSize(6), tBlockSize(6), color(_color) {
	this->initialize_index_lattice();
}

BlockDiracOperator::~BlockDiracOperator() { }

BlockDiracOperator* BlockDiracOperator::getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters, Color _color) {
	if (power == 1) {
		if (name == "DiracWilson") {
			BlockDiracWilsonOperator* result = new BlockDiracWilsonOperator(_color);
			result->setKappa(parameters.get<double>("kappa"));
			return result;
		} else if (name == "Improved") {
			BlockImprovedDiracWilsonOperator* result = new BlockImprovedDiracWilsonOperator(_color);
			result->setKappa(parameters.get<double>("kappa"));
			result->setCSW(parameters.get<double>("csw"));
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

void BlockDiracOperator::project(reduced_dirac_vector_t& output) {
#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		if (this->index_lattice[site] == 0) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_zero(output[site][mu]);
			}
		}
	}
}

void BlockDiracOperator::setBlockSize(const std::vector<unsigned int>& _blockSize) {
	xBlockSize = _blockSize[0];
	yBlockSize = _blockSize[1];
	zBlockSize = _blockSize[2];
	tBlockSize = _blockSize[3];
	this->initialize_index_lattice();
}

std::vector<unsigned int> BlockDiracOperator::getBlockSize() const {
	std::vector<unsigned int> result;
	result.push_back(xBlockSize);
	result.push_back(yBlockSize);
	result.push_back(zBlockSize);
	result.push_back(tBlockSize);
	return result;
}

void BlockDiracOperator::initialize_index_lattice() {
	typedef reduced_index_lattice_t::Layout Layout;
	int c;
	if (color == Black) c = 0;
	else c = 1;
#pragma omp parallel for
	for (int site = 0; site < index_lattice.localsize; ++site) {
		int is_even = (Layout::globalIndexX(site)/xBlockSize) + (Layout::globalIndexY(site)/yBlockSize)
			+ (Layout::globalIndexZ(site)/zBlockSize) + (Layout::globalIndexT(site)/tBlockSize);
		if ((is_even % 2) == c) this->index_lattice[site] = 1;
		else this->index_lattice[site] = 0;
	}
	index_lattice.updateHalo();
}

} /* namespace Update */
