/*
 * BlockDiracOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef BLOCKDIRACOPERATOR_H_
#define BLOCKDIRACOPERATOR_H_
#include "DiracOperator.h"
#include "../RandomSeed.h"

namespace Update {

class BlockDiracOperator : public DiracOperator {
public:
	BlockDiracOperator();
	BlockDiracOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa);
	~BlockDiracOperator();

	static DiracOperator* getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters);

	virtual void setBlockSize(const std::vector<unsigned int>& _blockSize);
	std::vector<unsigned int> getBlockSize() const;
protected:
	int xBlockSize;
	int yBlockSize;
	int zBlockSize;
	int tBlockSize;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
