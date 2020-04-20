#ifndef BASICSQUAREDIRACWILSONOPERATOR_H_
#define BASICSQUAREDIRACWILSONOPERATOR_H_
#include "DiracOperator.h"

namespace Update {

class BasicSquareDiracWilsonOperator : public Update::DiracOperator {
public:
	BasicSquareDiracWilsonOperator();
	BasicSquareDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., bool _gamma5 = true);
	~BasicSquareDiracWilsonOperator();

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
private:
	reduced_dirac_vector_t tmp;
};

} /* namespace Update */
#endif /* BASICSQUAREDIRACWILSONOPERATOR_H_ */
