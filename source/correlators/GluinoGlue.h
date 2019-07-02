#ifndef GLUINOGLUE_H
#define GLUINOGLUE_H

#include "LatticeSweep.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/BiConjugateGradient.h"
#include "fermion_measurements/StochasticEstimator.h"

namespace Update {

class GluinoGlue : public Update::LatticeSweep, public StochasticEstimator {
public:
	GluinoGlue();
	GluinoGlue(const GluinoGlue& toCopy);
	~GluinoGlue();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);

private:
	extended_dirac_vector_t source;
	extended_dirac_vector_t rho;
	extended_dirac_vector_t eta;
	extended_dirac_vector_t psi[4];
	extended_dirac_vector_t randomNoise;
	DiracOperator* diracOperator;
	BiConjugateGradient* biConjugateGradient;

	GaugeGroup cloverPlaquette(const extended_gauge_lattice_t& lattice, int site, int mu, int nu);
};

}

#endif
