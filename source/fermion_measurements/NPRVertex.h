#ifndef NPRVERTEX_H_
#define NPRVERTEX_H_

#include "LatticeSweep.h"
#include "fermion_measurements/StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/Solver.h"
#include "multigrid/MultiGridStochasticEstimator.h"
#include "utils/Gamma.h"

namespace Update {

class NPRVertex : public Update::LatticeSweep, public Update::MultiGridStochasticEstimator {
public:
	NPRVertex();
	NPRVertex(const NPRVertex& toCopy);
	~NPRVertex();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description&);

private:
	extended_dirac_vector_t randomNoise;
	extended_dirac_vector_t eta;
	extended_dirac_vector_t inverse_source[4*diracVectorLength];
	extended_dirac_vector_t source[4*diracVectorLength];
	DiracOperator* diracOperator;
	Solver* inverter;

	Gamma gamma;
};

} /* namespace Update */
#endif /* MESONCORRELATOR_H_ */
