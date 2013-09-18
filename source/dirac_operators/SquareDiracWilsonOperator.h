/*
 * SquareDiracWilsonOperator.h
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#ifndef SQUAREDIRACWILSONOPERATOR_H_
#define SQUAREDIRACWILSONOPERATOR_H_

#include "DiracOperator.h"
#include "DiracWilsonOperator.h"

namespace Update {

class SquareDiracWilsonOperator : public Update::DiracOperator {
public:
	SquareDiracWilsonOperator();
	SquareDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., bool _gamma5 = true);
	~SquareDiracWilsonOperator();

	/**
	 * This routine multiplies the DiracWilson operator two times to input and stores the result in output
	 * @param output
	 * @param input
	 */
	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	/**
	 * This routine multiplies the DiracWilson operator two times to vector1 and stores the result in output adding to it alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha);

	virtual FermionForce* getForce() const;

	virtual void setKappa(real_t _kappa);

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);
private:
	DiracWilsonOperator diracWilsonOperator;
	
	reduced_dirac_vector_t tmp;
};

} /* namespace Update */
#endif /* SQUAREDIRACWILSONOPERATOR_H_ */
