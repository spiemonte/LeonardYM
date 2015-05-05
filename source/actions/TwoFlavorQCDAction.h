/*
 * TwoFlavorAction.h
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#ifndef TWOFLAVORACTION_H_
#define TWOFLAVORACTION_H_

#include "GaugeForce.h"
#include "TwoFlavorFermionAction.h"
#include "GaugeAction.h"

namespace Update {

class TwoFlavorAction : public Energy, public Force {
public:
	TwoFlavorAction(GaugeAction* _gaugeAction, TwoFlavorFermionAction* _fermionAction);
	~TwoFlavorAction();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;

	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);

	virtual long_real_t energy(const environment_t& env);

	void setGaugeAction(GaugeAction* _gaugeAction);
	GaugeAction* getGaugeAction() const;

	void setFermionAction(TwoFlavorFermionAction* _fermionAction);
	TwoFlavorFermionAction* getFermionAction() const;
private:
	//The gluon part of the action
	GaugeAction* gaugeAction;
	//The fermion part of the action
	TwoFlavorFermionAction* fermionAction;
};

} /* namespace Update */
#endif /* TWOFLAVORQCDACTION_H_ */
