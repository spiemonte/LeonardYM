/*
 * NFlavorQCDAction.cpp
 *
 *  Created on: May 10, 2012
 *      Author: spiem_01
 */

#include "NFlavorQCDAction.h"
#include "inverters/MultishiftSolver.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

NFlavorAction::NFlavorAction(GaugeAction* _gaugeAction, NFlavorFermionAction* _fermionAction) : gaugeAction(_gaugeAction), fermionAction(_fermionAction) { }

NFlavorAction::~NFlavorAction() {
	delete gaugeAction;
	delete fermionAction;
}

GaugeGroup NFlavorAction::force(const environment_t& env, int site, int mu) const {
	//We return only the gauge force, fermion is already calculated
	return gaugeAction->force(env, site, mu);
}

void NFlavorAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	//Calculate the fermion force
	fermionAction->updateForce(forceLattice, env);

	//Add the gauge force
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			forceLattice[site][mu] += this->force(env, site, mu);
		}
	}

	forceLattice.updateHalo();//TODO is needed?
}

long_real_t NFlavorAction::energy(const environment_t& env) {
	return fermionAction->energy(env) + gaugeAction->energy(env);
}

void NFlavorAction::setGaugeAction(GaugeAction* _gaugeAction) {
	gaugeAction = _gaugeAction;
}

GaugeAction* NFlavorAction::getGaugeAction() const {
	return gaugeAction;
}

void NFlavorAction::setFermionAction(NFlavorFermionAction* _fermionAction) {
	fermionAction = _fermionAction;
}

NFlavorFermionAction* NFlavorAction::getFermionAction() const {
	return fermionAction;
}

} /* namespace Update */
