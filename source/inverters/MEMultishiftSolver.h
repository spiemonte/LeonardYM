#ifndef MEMULTISHIFTSOLVER_H_
#define MEMULTISHIFTSOLVER_H_
#include "MultishiftSolver.h"

namespace Update {

class MEMultishiftSolver : public MultishiftSolver {
public:
	MEMultishiftSolver(real_t _epsilon = 0.00000001, unsigned int _maxSteps = 3000);
	virtual ~MEMultishiftSolver();

	/**
	 * This function implements the multishift solver for the operator dirac (NB: it must hermitian and definite positive)
	 * @param dirac the dirac operator
	 * @param source
	 * @param solutions
	 * @param shifts
	 * @return false if the solver fails
	 */
	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, std::vector<extended_dirac_vector_t>& solutions, const std::vector<real_t>& shifts);

private:
	reduced_dirac_vector_t residual;
	reduced_dirac_vector_t p;
	reduced_dirac_vector_t tmp;
};

} /* namespace Update */
#endif /* MEMULTISHIFTSOLVER_H_ */
