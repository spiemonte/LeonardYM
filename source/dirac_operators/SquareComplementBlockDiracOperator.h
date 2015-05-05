/*
 * SquareComplementBlockDiracOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef SQUARECOMPLEMENTBLOCKDIRACOPERATOR_H_
#define SQUARECOMPLEMENTBLOCKDIRACOPERATOR_H_
#include "DiracWilsonOperator.h"
#include "ComplementBlockDiracWilsonOperator.h"
#include "SquareBlockDiracWilsonOperator.h"
#include "../BiConjugateGradient.h"

namespace Update {

class SquareComplementBlockDiracOperator : public BlockDiracOperator {
public:
	SquareComplementBlockDiracOperator(ComplementBlockDiracOperator* _K);
	~SquareComplementBlockDiracOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual void setKappa(real_t _kappa);

	//Set the precision of the inner inverter
	void setPrecision(const real_t& _precision);
	void setMaximumSteps(int steps);
private:
	ComplementBlockDiracOperator* K;
	reduced_dirac_vector_t tmpVector;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
