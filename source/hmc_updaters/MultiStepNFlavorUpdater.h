#ifndef MULTISTEPNFLAVORUPDATER_H_
#define MULTISTEPNFLAVORUPDATER_H_

#include "LatticeSweep.h"
#include "FermionHMCUpdater.h"
#include "dirac_functions/RationalApproximation.h"
#include "dirac_operators/BlockDiracOperator.h"
#include "actions/NFlavorQCDAction.h"
#include "actions/GaugeAction.h"

#include <vector>

namespace Update {

class MultiStepNFlavorUpdater: public Update::LatticeSweep, public Update::FermionHMCUpdater {
public:
	MultiStepNFlavorUpdater();
	MultiStepNFlavorUpdater(const MultiStepNFlavorUpdater& toCopy);
	~MultiStepNFlavorUpdater();

	virtual void execute(environment_t& environment);

	static void registerParameters(po::options_description& desc);

protected:
	void initializeApproximations(environment_t& environment);

	void checkTheory(const environment_t& environment) const;

protected:
	//The new environment, provided by HMC
	environment_t environmentNew;
	//The conjugate momenta
	extended_gauge_lattice_t momenta;
	//The pseudofermion field
	std::vector<extended_dirac_vector_t> pseudofermions;
	std::vector<RationalApproximation> rationalApproximationsHeatBath;
	std::vector<RationalApproximation> rationalApproximationsMetropolis;
	std::vector< std::vector<RationalApproximation> > rationalApproximationsForce;
	//The vector of polynomia
	//The tmp pseudofermion field
	extended_dirac_vector_t tmp_pseudofermion;
	//The action of the theory
	NFlavorAction* nFlavorAction;
	//The gauge action of the theory
	GaugeAction* gaugeAction;
	//The fermion action of the theory
	NFlavorFermionAction** fermionAction;
	//The dirac operators for the theory
	DiracOperator* squareDiracOperatorMetropolis;
	DiracOperator* diracOperatorMetropolis;
	DiracOperator* squareDiracOperatorForce;
	DiracOperator* diracOperatorForce;

	MultishiftSolver* multishiftSolver;
	BlockDiracOperator* blackBlockDiracOperator;
	BlockDiracOperator* redBlockDiracOperator;
};

} /* namespace Update */
#endif /* NFLAVORUPDATER_H_ */
