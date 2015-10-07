#ifndef SINGLETOPERATORS_H
#define SINGLETOPERATORS_H
#include "StochasticEstimator.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/BiConjugateGradient.h"

namespace Update {

class SingletOperators : public StochasticEstimator, LatticeSweep, WilsonFlow {
public:
	SingletOperators();

	virtual void execute(environment_t& environment);

private:
	extended_dirac_vector_t randomNoise;
	extended_dirac_vector_t tmp_square;
	extended_dirac_vector_t tmp;
	DiracOperator* squareDiracOperator;
	DiracOperator* diracOperator;
	BiConjugateGradient* biConjugateGradient;
};

}

#endif
