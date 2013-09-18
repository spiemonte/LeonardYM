/*
 * BasicDiracWilsonOperator.h
 *
 *  Created on: Feb 25, 2013
 *      Author: spiem_01
 */

#ifndef BASICDIRACWILSONOPERATOR_H_
#define BASICDIRACWILSONOPERATOR_H_
#include "DiracOperator.h"

namespace Update {

class BasicDiracWilsonOperator : public DiracOperator {
public:
	BasicDiracWilsonOperator();
	BasicDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., bool _gamma5 = true);
	virtual ~BasicDiracWilsonOperator();

	/**
	 * This routine multiplies the DiracWilson operator to input and stores the result in output
	 * @param output
	 * @param input
	 */
	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	/**
	 * This routine multiplies the DiracWilson operator to vector1 and stores the result in output adding to it alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha);

	virtual FermionForce* getForce() const;
private:
	BasicDiracWilsonOperator(const DiracOperator&) { }
};

} /* namespace Update */
#endif /* BASICDIRACWILSONOPERATOR_H_ */
