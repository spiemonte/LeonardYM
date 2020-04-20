#ifndef SQUAREIMPROVEDDIRACWILSONOPERATOR_H_
#define SQUAREIMPROVEDDIRACWILSONOPERATOR_H_
#include "DiracOperator.h"
#include "ImprovedDiracWilsonOperator.h"

namespace Update {

class SquareImprovedDiracWilsonOperator : public Update::DiracOperator  {
public:
	SquareImprovedDiracWilsonOperator();
	SquareImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., real_t _csw = 1., bool _gamma5 = true);
	~SquareImprovedDiracWilsonOperator();

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

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual FermionForce* getForce() const;

	real_t getCSW() const;
	void setCSW(real_t _csw);

	virtual void setKappa(real_t _kappa);
private:
	ImprovedDiracWilsonOperator improvedDiracWilsonOperator;
	reduced_dirac_vector_t tmp;
	real_t csw;
};

} /* namespace Update */
#endif /* SQUAREIMPROVEDDIRACWILSONOPERATOR_H_ */
