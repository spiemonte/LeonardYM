#ifndef OVERLAPOPERATOR_H_
#define OVERLAPOPERATOR_H_
#include "DiracOperator.h"
#include "DiracWilsonOperator.h"
#include "SquareDiracWilsonOperator.h"
#include "dirac_functions/Polynomial.h"

namespace Update {

class OverlapOperator : public DiracOperator {
public:
	OverlapOperator();
	OverlapOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., bool _gamma5 = true);
	virtual ~OverlapOperator();

	/**
	 * This routine multiplies the Overlap operator to input and stores the result in output
	 * @param output
	 * @param input
	 */
	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	/**
	 * This routine multiplies the Overlap operator to vector1 and stores the result in output adding to it alpha*vector2
	 * @param output
	 * @param vector1
	 * @param vector2
	 * @param alpha
	 */
	virtual void multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha);

	virtual void setKappa(real_t _kappa);
	void setMass(real_t _mass);
	real_t getMass() const;

	void setSquareRootApproximation(const Polynomial& _squareRootApproximation);
	Polynomial& getSquareRootApproximation();
	const Polynomial& getSquareRootApproximation() const;

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	virtual FermionForce* getForce() const;
protected:
	real_t mass;

	Polynomial squareRootApproximation;
	DiracWilsonOperator diracWilsonOperator;
	SquareDiracWilsonOperator squareDiracWilsonOperator;

	reduced_dirac_vector_t tmp1;
	reduced_dirac_vector_t tmp2;
private:
	OverlapOperator(const DiracOperator&) { }
};

} /* namespace Update */
#endif /* OVERLAPOPERATOR_H_ */
