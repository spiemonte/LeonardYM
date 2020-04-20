#ifndef CHRONOLOGICALMULTISHIFTSOLVER_H_
#define CHRONOLOGICALMULTISHIFTSOLVER_H_
#include "MultishiftSolver.h"

namespace Update {

class ChronologicalMultishiftSolver : public MultishiftSolver {
public:
	ChronologicalMultishiftSolver(double _epsilon = 0.00000001, unsigned int _maxSteps = 3000);
	~ChronologicalMultishiftSolver();

	/**
	 * This function implements the multishift solver for the operator dirac (NB: it must hermitian and definite positive)
	 * The algorithm works in a chronological way, supposing that solutions is the vector that is set to the old solutions
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

	reduced_dirac_vector_t initial_guess;
};

} /* namespace Update */
#endif /* CHRONOLOGICALMULTISHIFTSOLVER_H_ */
