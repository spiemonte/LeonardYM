/*
 * BandTwoFlavorGlueballMeasure.h
 *
 *  Created on: Jan 21, 2013
 *      Author: spiem_01
 */

#ifndef BANDTWOFLAVORUPDATER_H_
#define BANDTWOFLAVORUPDATER_H_

#include "LatticeSweep.h"
#include "FermionHMCUpdater.h"
#include "BandAction.h"
#include "TwoFlavorQCDAction.h"

namespace Update {

class BandTwoFlavorUpdater : public Update::LatticeSweep, public Update::FermionHMCUpdater {
public:
	BandTwoFlavorUpdater();
	~BandTwoFlavorUpdater();

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
	//The band action used, it updates only the links in the bands
	BandAction* bandAction;
	BandAction* bandGaugeAction;
	BandAction* bandFermionAction;
};

} /* namespace Update */
#endif /* BANDTWOFLAVORGLUEBALLMEASURE_H_ */
