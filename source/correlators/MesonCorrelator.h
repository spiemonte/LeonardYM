#ifndef MESONCORRELATOR_H_
#define MESONCORRELATOR_H_
#include "LatticeSweep.h"
#include "fermion_measurements/StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/Solver.h"
#include <map>

namespace Update {

class MesonCorrelator: public Update::LatticeSweep, public Update::StochasticEstimator {
public:
	MesonCorrelator();
	MesonCorrelator(const MesonCorrelator& toCopy);
	~MesonCorrelator();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);

private:
	reduced_dirac_vector_t randomNoise;
	reduced_dirac_vector_t propagator[4*diracVectorLength];
	reduced_dirac_vector_t tmp;
	DiracOperator* diracOperator;
	Solver* inverter;
};

} /* namespace Update */
#endif /* MESONCORRELATOR_H_ */
