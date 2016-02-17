#ifndef MULTiGRIDCONJUGATEGRADIENT_H
#define MULTIGRIDCONJUGATEGRADIENT_H
#include "MultiGridVector.h"
#include "MultiGridOperator.h"

namespace Update {

class MultiGridConjugateGradientSolver {
public:
	MultiGridConjugateGradientSolver();

	bool solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source_hat, multigrid_vector_t& solution_hat);

	void setPrecision(const real_t& _precision);

	void setMaximumSteps(int _maxSteps);
private:
	real_t precision;
	int maxSteps;

	multigrid_vector_t r_hat;
	multigrid_vector_t p_hat;
	multigrid_vector_t tmp_hat;
};

}

#endif

