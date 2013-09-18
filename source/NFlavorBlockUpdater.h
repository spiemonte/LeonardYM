/*
 * NFlavorQCDUpdater.h
 *
 *  Created on: May 2, 2012
 *      Author: spiem_01
 */

#ifndef NFLAVORBLOCKUPDATER_H_
#define NFLAVORBLOCKUPDATER_H_

#include "LatticeSweep.h"
#include "FermionHMCUpdater.h"
#include "Polynomial.h"
#include "RationalApproximation.h"
#include "NFlavorQCDUpdater.h"
#include "GaugeAction.h"

#include <vector>

namespace Update {

class NFlavorBlockUpdater: public Update::NFlavorQCDUpdater {
public:
	NFlavorBlockUpdater();
	NFlavorBlockUpdater(const NFlavorBlockUpdater& toCopy);
	~NFlavorBlockUpdater();

	virtual void execute(environment_t& environment);
protected:
	void initializeCorrectionStepApproximations(const environment_t& environment);
	void checkCorrectionStepApproximations(const environment_t& environment) const;

	void stochasticCorrectionStep(environment_t& environment, const environment_t& environmentNew);
	
private:
	std::vector<Polynomial> polynomialApproximationsInverse;
	std::vector<Polynomial> polynomialApproximationsDirect;
};

} /* namespace Update */
#endif /* NFLAVORQCDUPDATER_H_ */
