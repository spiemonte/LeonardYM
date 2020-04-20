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

#ifdef ENABLE_MPI
	void evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& output, const extended_dirac_vector_t& input);
#endif
	void evaluate(DiracOperator* diracOperator, reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input);	

	void evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& output, const extended_dirac_vector_t& input, DiracOperator* preconditioner);

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
