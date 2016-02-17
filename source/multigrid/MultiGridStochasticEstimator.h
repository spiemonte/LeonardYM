#ifndef MULTIGRIDSTOCHASTICESTIMATOR_H
#define MULTIGRIDSTOCHASTICESTIMATOR_H
#include "BlockBasis.h"
#include "MultiGridBiConjugateGradient.h"
#include "inverters/Solver.h"
#include "fermion_measurements/StochasticEstimator.h"

namespace Update {

class MultiGridStochasticEstimator : public StochasticEstimator {
public:
	MultiGridStochasticEstimator();
	MultiGridStochasticEstimator(BlockBasis* _blockBasis);
	MultiGridStochasticEstimator(const MultiGridStochasticEstimator& toCopy);
	~MultiGridStochasticEstimator();

	void setBlockBasis(BlockBasis* _blockBasis);

	void getMultigridVectors(DiracOperator* dirac, extended_dirac_vector_t& mg_source, extended_dirac_vector_t& inverse_mg_source);

	void getResidualVectors(Solver* solver, DiracOperator* dirac, extended_dirac_vector_t& residual_source, extended_dirac_vector_t& inverse_residual_source);
private:
	BlockBasis* blockBasis;

	MultiGridBiConjugateGradientSolver* biMgSolver;
};

}


#endif

