#ifndef XSPACECORRELATORS_H
#define XSPACECORRELATORS_H
#include "StochasticEstimator.h"
#include "LatticeSweep.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/BiConjugateGradient.h"
#include "wilson_flow/WilsonFlow.h"

namespace Update {

class XSpaceCorrelators : public StochasticEstimator, public WilsonFlow {
public:
	XSpaceCorrelators();

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
