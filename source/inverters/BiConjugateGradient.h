/*
 * BiConjugateGradient.h
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#ifndef BICONJUGATEGRADIENT_H_
#define BICONJUGATEGRADIENT_H_
#include "dirac_operators/DiracOperator.h"
#include "Solver.h"

namespace Update {

class BiConjugateGradient : public Solver {
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
	using Solver::solve;

	BiConjugateGradient();
	~BiConjugateGradient();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess = 0);
	virtual bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0);
	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, const std::complex<real_t>& alpha, reduced_dirac_vector_t const* initial_guess = 0);
	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, int l, reduced_dirac_vector_t const* initial_guess = 0);

};

} /* namespace Update */
#endif /* BICONJUGATEGRADIENT_H_ */
