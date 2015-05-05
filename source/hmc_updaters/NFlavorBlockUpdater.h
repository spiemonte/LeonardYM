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
#include "dirac_functions/RationalApproximation.h"
#include "MultiStepNFlavorQCDUpdater.h"
#include "actions/GaugeAction.h"

#include <vector>

namespace Update {

class BiConjugateGradient;

class NFlavorBlockUpdater: public Update::MultiStepNFlavorQCDUpdater {
public:
	NFlavorBlockUpdater();
	NFlavorBlockUpdater(const NFlavorBlockUpdater& toCopy);
	~NFlavorBlockUpdater();

	virtual void execute(environment_t& environment);
protected:
	long_real_t logDeterminant(const environment_t& environment, real_t alpha);
	long_real_t logDeterminant_test(const environment_t& environment, real_t alpha);
	long_real_t logDeterminant(const environment_t& environment_new, const environment_t& environment_old, real_t alpha);

	reduced_dirac_vector_t temp1, temp2, temp3, temp4, result;
	long_real_t log_determinant;
	BiConjugateGradient* biConjugateGradient;
};

} /* namespace Update */
#endif /* NFLAVORQCDUPDATER_H_ */
