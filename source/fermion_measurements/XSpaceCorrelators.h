#ifndef XSPACECORRELATORS_H
#define XSPACECORRELATORS_H
#include "StochasticEstimator.h"
#include "LatticeSweep.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/BiConjugateGradient.h"
#include "wilson_flow/WilsonFlow.h"
#include "dirac_operators/GammaOperators.h"

namespace Update {

class XSpaceCorrelators : public StochasticEstimator, public WilsonFlow {
public:
	XSpaceCorrelators();

	virtual void execute(environment_t& environment);

private:
	extended_dirac_vector_t tmp, source;
	DiracOperator* squareDiracOperator;
	DiracOperator* diracOperator;
	BiConjugateGradient* biConjugateGradient;
	
	GammaOperators gammaOperators;
	Gamma gammas;
};

}

#endif
