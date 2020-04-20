#ifndef SINGLETOPERATORS_H_
#define SINGLETOPERATORS_H_

#include "LatticeSweep.h"
#include "fermion_measurements/StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/Solver.h"
#include "multigrid/MultiGridStochasticEstimator.h"
#include "utils/Gamma.h"

namespace Update {

class SingletOperators : public Update::LatticeSweep, public Update::MultiGridStochasticEstimator {
public:
	SingletOperators();
	SingletOperators(const SingletOperators& toCopy);
	~SingletOperators();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description&);

private:
	extended_dirac_vector_t randomNoise;
	extended_dirac_vector_t eta;
	DiracOperator* diracOperator;
	DiracOperator* squareDiracOperator;
	Solver* inverter;

	Gamma gamma;
};

} /* namespace Update */
#endif /* MESONCORRELATOR_H_ */
