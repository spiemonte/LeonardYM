/*
 * MesonCorrelator.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef MESONCORRELATOR_H_
#define MESONCORRELATOR_H_

#include "LatticeSweep.h"
#include "fermion_measurements/StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/Solver.h"

namespace Update {

class MesonCorrelator: public Update::LatticeSweep, public Update::StochasticEstimator {
public:
	MesonCorrelator();
	MesonCorrelator(const MesonCorrelator& toCopy);
	~MesonCorrelator();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);

private:
	extended_dirac_vector_t randomNoise;
	extended_dirac_vector_t tmp[4*diracVectorLength];
	DiracOperator* diracOperator;
	Solver* inverter;
};

} /* namespace Update */
#endif /* MESONCORRELATOR_H_ */
