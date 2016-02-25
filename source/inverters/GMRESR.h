/*
 * GMRESR.h
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#ifndef GMRESR_H_
#define GMRESR_H_
#include "dirac_operators/DiracOperator.h"
#include "Solver.h"

namespace Update {

class GMRESR : public Solver {
	unsigned int conjugateSpaceDimension;

	reduced_dirac_vector_t* c;
	reduced_dirac_vector_t* u;
	reduced_dirac_vector_t r;
	
public:
	using Solver::solve;

	GMRESR();
	GMRESR(const GMRESR& toCopy);
	~GMRESR();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess = 0);

	virtual bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0) {
		return this->solve(dirac, source, solution, 0, initial_guess);
	}

#ifdef ENABLE_MPI
	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, DiracOperator* preconditioner, extended_dirac_vector_t const* initial_guess = 0);
#endif

	void setConjugateSpaceDimension(unsigned int _maxSteps);
	unsigned int getConjugateSpaceDimension() const;
};

} /* namespace Update */
#endif /* BICONJUGATEGRADIENT_H_ */
