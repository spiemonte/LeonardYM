/*
 * PureGaugeHMCUpdater.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef PUREGAUGEHMCUPDATER_H_
#define PUREGAUGEHMCUPDATER_H_

#include "LatticeSweep.h"
#include "Integrate.h"
#include "HMCUpdater.h"

namespace Update {

class PureGaugeHMCUpdater: public Update::LatticeSweep, public HMCUpdater {
public:
	PureGaugeHMCUpdater();
	~PureGaugeHMCUpdater();

	virtual void execute(environment_t& environment);

private:
	//The new environment, provided by HMC
	environment_t environmentNew;
	//The conjugate momenta
	extended_gauge_lattice_t momenta;
};

} /* namespace Update */
#endif /* PUREGAUGEHMCUPDATER_H_ */
