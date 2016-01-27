/*
 * SquareComplementBlockDiracWilsonOperator.h
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#ifndef SQUARECOMPLEMENTBLOCKDIRACWILSONOPERATOR_H_
#define SQUARECOMPLEMENTBLOCKDIRACWILSONOPERATOR_H_
#include "DiracWilsonOperator.h"
#include "ComplementBlockDiracOperator.h"
#include "SquareBlockDiracWilsonOperator.h"
#include "inverters/BiConjugateGradient.h"

namespace Update {

class SquareComplementBlockDiracWilsonOperator : public DiracOperator {
public:
	SquareComplementBlockDiracWilsonOperator();
	SquareComplementBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~SquareComplementBlockDiracWilsonOperator();

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha);

	virtual FermionForce* getForce() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual void setKappa(real_t _kappa);
	
	void setBlockSize(const std::vector<unsigned int>&_blockSize);
	std::vector<unsigned int> getBlockSize() const;

	//Set the precision of the inner inverter
	void setPrecision(const real_t& _precision);
	real_t getPrecision() const;
private:
	DiracWilsonOperator diracWilsonOperator;
	BlockDiracWilsonOperator blockDiracWilsonOperator;
	SquareBlockDiracWilsonOperator squareBlockDiracWilsonOperator;
	BiConjugateGradient biConjugateGradient;
	reduced_dirac_vector_t tmpVector;
	bool log;
	int counterSteps;
};

} /* namespace Update */
#endif /* RANDOMDIRACWILSONOPERATOR_H_ */
