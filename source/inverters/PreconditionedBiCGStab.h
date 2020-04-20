#ifndef PRECONDITIONEDBICGSTAB_H_
#define PRECONDITIONEDBICGSTAB_H_
#include "dirac_operators/DiracOperator.h"
#include "BiConjugateGradient.h"
#include "Solver.h"

namespace Update {

class PreconditionedBiCGStab : public Solver {
public:
	using Solver::solve;

	PreconditionedBiCGStab(const PreconditionedBiCGStab& toCopy);
	PreconditionedBiCGStab();
	~PreconditionedBiCGStab();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess = 0);

	virtual bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0) {
		return this->solve(dirac, source, solution, 0, initial_guess);
	}

#ifdef ENABLE_MPI
	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, DiracOperator* preconditioner, extended_dirac_vector_t const* initial_guess = 0);
#endif
	
	void setUseEvenOddPreconditioning(bool);
private:
	bool useEvenOddPreconditioning;
	BiConjugateGradient* biConjugateGradient;
};

} /* namespace Update */
#endif /* BICONJUGATEGRADIENT_H_ */
