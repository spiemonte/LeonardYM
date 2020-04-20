#ifndef MULTIGRIDMEMULTISHIFTSOLVER_H_
#define MULTIGRIDMEMULTISHIFTSOLVER_H_
#include "MultishiftSolver.h"
#include "multigrid/MultiGridSolver.h"
#include "dirac_operators/TwistedDiracOperator.h"

namespace Update {

class MultiGridMEMultishiftSolver : public MultishiftSolver, public MultiGridSolver {
public:
	MultiGridMEMultishiftSolver(int basisDimension, const std::vector<unsigned int>& _blockSize, BlockDiracOperator* _blackBlockDiracOperator, BlockDiracOperator* _redBlockDiracOperator, real_t _epsilon = 0.00000001, unsigned int _maxSteps = 3000);
	MultiGridMEMultishiftSolver(const MultiGridMEMultishiftSolver& toCopy);
	virtual ~MultiGridMEMultishiftSolver();

	/**
	 * This function implements the multishift solver for the operator dirac (NB: it must hermitian and definite positive)
	 * @param dirac the dirac operator
	 * @param source
	 * @param solutions
	 * @param shifts
	 * @return false if the solver fails
	 */
	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, std::vector<extended_dirac_vector_t>& solutions, const std::vector<real_t>& shifts);
	
	using MultiGridSolver::solve;
private:
	reduced_dirac_vector_t p;
	reduced_dirac_vector_t r;
};

} /* namespace Update */
#endif /* MEMULTISHIFTSOLVER_H_ */
