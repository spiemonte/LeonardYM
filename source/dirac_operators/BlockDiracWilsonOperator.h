/*
 * BlockDiracWilsonOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef BLOCKDIRACWILSONOPERATOR_H_
#define BLOCKDIRACWILSONOPERATOR_H_
#include "BlockDiracOperator.h"
#include "../RandomSeed.h"

namespace Update {

class BlockDiracWilsonOperator : public BlockDiracOperator {
public:
	BlockDiracWilsonOperator();
	BlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~BlockDiracWilsonOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
