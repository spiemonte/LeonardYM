/*
 * ComplementBlockDiracWilsonOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef COMPLEMENTBLOCKDIRACWILSONOPERATOR_H_
#define COMPLEMENTBLOCKDIRACWILSONOPERATOR_H_
#include "DiracWilsonOperator.h"
#include "BlockDiracOperator.h"

namespace Update {

class ComplementBlockDiracWilsonOperator : public BlockDiracOperator {
public:
	ComplementBlockDiracWilsonOperator();
	ComplementBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~ComplementBlockDiracWilsonOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual void setKappa(real_t _kappa);

	void setDirection(unsigned int mu);
	unsigned int getDirection() const;

private:
	DiracWilsonOperator diracWilsonOperator;
	unsigned int dir;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
