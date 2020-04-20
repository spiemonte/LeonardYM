#ifndef SQUAREBLOCKDIRACWILSONOPERATOR_H_
#define SQUAREBLOCKDIRACWILSONOPERATOR_H_

#include "DiracOperator.h"
#include "BlockDiracWilsonOperator.h"

namespace Update {

class SquareBlockDiracWilsonOperator : public BlockDiracOperator {
public:
	SquareBlockDiracWilsonOperator();
	SquareBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0.);
	~SquareBlockDiracWilsonOperator();

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

	virtual void setBlockSize(const std::vector<unsigned int>& _blockSize);
private:
	BlockDiracWilsonOperator blockDiracWilsonOperator;
	
	reduced_dirac_vector_t tmp;
};

} /* namespace Update */
#endif /* SQUAREDIRACWILSONOPERATOR_H_ */
