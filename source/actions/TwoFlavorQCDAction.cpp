#include "TwoFlavorQCDAction.h"

namespace Update {

TwoFlavorAction::TwoFlavorAction(GaugeAction* _gaugeAction, TwoFlavorFermionAction* _fermionAction) : gaugeAction(_gaugeAction), fermionAction(_fermionAction) { }

TwoFlavorAction::~TwoFlavorAction() {
	delete gaugeAction;
	delete fermionAction;
}

GaugeGroup TwoFlavorAction::force(const environment_t& env, int site, int mu) const {
	//We return back only the gauge force, the fermion part is already computed
	return gaugeAction->force(env, site, mu);
}

void TwoFlavorAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	//First calculate the fermion force
	fermionAction->updateForce(forceLattice, env);
	//Then add the gauge force
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			forceLattice[site][mu] += this->force(env, site, mu);
		}
	}
	forceLattice.updateHalo();//TODO is needed?
}

long_real_t TwoFlavorAction::energy(const environment_t& env) {
	return fermionAction->energy(env)+gaugeAction->energy(env);
}

void TwoFlavorAction::setGaugeAction(GaugeAction* _gaugeAction) {
	gaugeAction = _gaugeAction;
}

GaugeAction* TwoFlavorAction::getGaugeAction() const {
	return gaugeAction;
}

void TwoFlavorAction::setFermionAction(TwoFlavorFermionAction* _fermionAction) {
	fermionAction = _fermionAction;
}

TwoFlavorFermionAction* TwoFlavorAction::getFermionAction() const {
	return fermionAction;
}

} /* namespace Update */
