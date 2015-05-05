/*
 * RationalApproximation.h
 *
 *  Created on: May 10, 2012
 *      Author: spiem_01
 */

#ifndef RATIONALAPPROXIMATION_H_
#define RATIONALAPPROXIMATION_H_
#include "Environment.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/MultishiftSolver.h"
#include "inverters/BiConjugateGradient.h"
#include <vector>

namespace Update {

class RationalApproximation {
public:
	RationalApproximation();
	RationalApproximation(const RationalApproximation& toCopy);
	RationalApproximation(const std::vector< real_t >& _alphas, const std::vector< real_t >& _betas);
	~RationalApproximation();

	void evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& output, const extended_dirac_vector_t& input, DiracOperator* preconditioner = 0);

	complex evaluate(const complex& x) const;

	void setAlphas(const std::vector< real_t >& _alphas);
	void setBetas(const std::vector< real_t >& _betas);
	std::vector< real_t > getAlphas() const;
	std::vector< real_t > getBetas() const;

	void setPrecision(const real_t& _precision);
	real_t getPrecision() const;

	void setMaximumRecursion(unsigned int maximumRecursion);
	unsigned int getMaximumRecursion() const;

	void setPreconditionerRecursion(unsigned int recursion);
	unsigned int getPreconditionerRecursion() const;

	void setPreconditionerPrecision(const real_t& precision);
	real_t getPreconditionerPrecision() const;

private:
	std::vector< real_t > alphas;
	std::vector< real_t > betas;
	//This is the tmp memory used for the evaluation of the RationalApproximation
	static std::vector< extended_dirac_vector_t > tmp;
	extended_dirac_vector_t tmp1;
	//The precision of the inverter
	real_t precision;
	//The precision of the preconditioner
	real_t preconditionerPrecision;
	//The maximum steps for the inverter
	unsigned int maximumSteps;


	MultishiftSolver* multishiftSolver;
	BiConjugateGradient* biConjugateGradient;

	unsigned int preconditionerRecursion;
};

} /* namespace Update */
#endif /* RATIONALAPPROXIMATION_H_ */
