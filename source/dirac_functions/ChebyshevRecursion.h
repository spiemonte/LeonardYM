/*
 * Polynomial.h
 *
 *  Created on: Apr 30, 2012
 *      Author: spiem_01
 */

#ifndef CHEBYSHEVRECURSION_H_
#define CHEBYSHEVRECURSION_H_
#include "Environment.h"
#include "dirac_operators/DiracOperator.h"
#include <vector>

namespace Update {

class ChebyshevRecursion {
public:
	ChebyshevRecursion();
	ChebyshevRecursion(real_t a, real_t b, unsigned int n);
	ChebyshevRecursion(const ChebyshevRecursion& copy);
	~ChebyshevRecursion();

	void evaluate(DiracOperator* diracOperator, reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);

	complex evaluate(const complex& x) const;

	void setParameters(real_t a, real_t b, unsigned int n);
private:
	real_t a;
	real_t b;
	unsigned int n;
};

} /* namespace Update */
#endif /* POLYNOMIAL_H_ */
