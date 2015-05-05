/*
 * BiConjugateGradient.h
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#ifndef BICONJUGATEGRADIENT_H_
#define BICONJUGATEGRADIENT_H_
#include "dirac_operators/DiracOperator.h"
#include <vector>

namespace Update {

class BiConjugateGradient {
	double epsilon;
	double lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;

	reduced_dirac_vector_t residual;
	reduced_dirac_vector_t residual_hat;
	reduced_dirac_vector_t p;
	reduced_dirac_vector_t nu;
	reduced_dirac_vector_t s;
	reduced_dirac_vector_t t;
	
	//Vectors needed for preconditioning
	reduced_dirac_vector_t p_tilde;
	reduced_dirac_vector_t s_tilde;
public:
	BiConjugateGradient();
	~BiConjugateGradient();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess = 0);
	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0);
	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, const std::complex<real_t>& alpha, reduced_dirac_vector_t const* initial_guess = 0);
	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, int l, reduced_dirac_vector_t const* initial_guess = 0);

#ifdef ENABLE_MPI
	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, extended_dirac_vector_t const* initial_guess = 0);
#endif

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;
};

} /* namespace Update */
#endif /* BICONJUGATEGRADIENT_H_ */
