#ifndef MULTIGRIDBICONJUGATEGRADIENT_H
#define MULTIGRIDBICONJUGATEGRADIENT_H
#include "MultiGridVector.h"
#include "MultiGridOperator.h"

namespace Update {

class MultiGridBiConjugateGradientSolver {
public:
	MultiGridBiConjugateGradientSolver();

	bool solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source, multigrid_vector_t& solution);

	void setPrecision(const real_t& _precision);

	void setMaximumSteps(unsigned int _maxSteps);
private:
	real_t precision;
	unsigned int maxSteps;

	multigrid_vector_t residual;
	multigrid_vector_t residual_hat;
	multigrid_vector_t p;
	multigrid_vector_t nu;
	multigrid_vector_t s;
	multigrid_vector_t t;
};

}

#endif
