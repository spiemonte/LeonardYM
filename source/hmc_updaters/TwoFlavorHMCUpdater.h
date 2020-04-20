#ifndef TWOFLAVORHMCUPDATER_H_
#define TWOFLAVORHMCUPDATER_H_

#include "FermionHMCUpdater.h"
#include "LatticeSweep.h"
#include "actions/TwoFlavorQCDAction.h"
#include "inverters/BiConjugateGradient.h"

namespace Update {

class TwoFlavorHMCUpdater: public Update::LatticeSweep, public Update::FermionHMCUpdater {
public:
	TwoFlavorHMCUpdater();
	TwoFlavorHMCUpdater(const TwoFlavorHMCUpdater& toCopy);
	~TwoFlavorHMCUpdater();

	virtual void execute(environment_t& environment);
private:
	//The new environment, provided by HMC
	environment_t environmentNew;
	//The conjugate momenta
	extended_gauge_lattice_t momenta;
	//The pseudofermion field
	extended_dirac_vector_t pseudofermion;
	//The tmp pseudofermion field
	extended_dirac_vector_t tmp_pseudofermion;
	//The action for the theory
	TwoFlavorAction* action;
	//the Gauge action for the theory
	GaugeAction* gaugeAction;
	//the fermion action for the theory
	TwoFlavorFermionAction* fermionAction;
	//The dirac operator for the theory
	DiracOperator* diracOperator;
	//The inverter
	BiConjugateGradient* solver;
};

} /* namespace Update */
#endif /* TWOFLAVORHCMUPDATER_H_ */
