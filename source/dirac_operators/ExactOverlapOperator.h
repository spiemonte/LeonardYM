/*
 * ExactOverlapOperator.h
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#ifndef EXACTOVERLAPOPERATOR_H_
#define EXACTOVERLAPOPERATOR_H_
#include "OverlapOperator.h"
#include "fermion_measurements/DiracEigenSolver.h"


namespace Update {

class ExactOverlapOperator : public OverlapOperator {
public:
	ExactOverlapOperator();
	ExactOverlapOperator(const ExactOverlapOperator& copy);
	ExactOverlapOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa = 0., bool _gamma5 = true);
	virtual ~ExactOverlapOperator();

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

	virtual void setLattice(const extended_fermion_lattice_t& _lattice);

	void setNumberOfEigenvalues(unsigned int _numberOfEigenvalues);
	unsigned int getNumberOfEigenvalues() const;

	DiracEigenSolver* getDiracEigenSolver() const;
protected:
	void computeEigenvalues();
	bool checkEigenvalues();

private:
	ExactOverlapOperator(const DiracOperator&) { }

	DiracEigenSolver* diracEigenSolver;
	bool recomputeEigenvalues;

	std::vector< real_t > computed_eigenvalues;
	std::vector< reduced_dirac_vector_t > computed_eigenvectors;

	unsigned int numberOfEigenvalues;
};

} /* namespace Update */
#endif /* OVERLAPOPERATOR_H_ */
