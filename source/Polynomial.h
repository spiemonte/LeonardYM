/*
 * Polynomial.h
 *
 *  Created on: Apr 30, 2012
 *      Author: spiem_01
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_
#include "Environment.h"
#include "dirac_operators/DiracOperator.h"
#include <vector>

namespace Update {

class Polynomial {
public:
	Polynomial();
	Polynomial(const std::vector< complex >& roots, const complex& scaling);
	~Polynomial();

	void evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& output, const extended_dirac_vector_t& input);

	complex evaluate(const complex& x) const;

	void setScaling(const complex& _scaling);
	complex getScaling() const;

	void setRoots(const std::vector<complex>& _roots);
	std::vector<complex> getRoots() const;

private:
	std::vector< complex > roots;
	complex scaling;
	//Static tmp pointer to the vector needed for the calculation of the polynomial
	reduced_dirac_vector_t tmp1;
	reduced_dirac_vector_t tmp2;
};

} /* namespace Update */
#endif /* POLYNOMIAL_H_ */
