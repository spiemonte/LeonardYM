#ifndef NFLAVORQCDACTION_H_
#define NFLAVORQCDACTION_H_
#include "hmc_forces/Force.h"
#include "GaugeAction.h"
#include "NFlavorFermionAction.h"

#include <vector>

namespace Update {

class NFlavorAction : public Energy, public Force {
public:
	NFlavorAction(GaugeAction* _gaugeAction, NFlavorFermionAction* _fermionAction);
	virtual ~NFlavorAction();

	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;

	virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);

	virtual long_real_t energy(const environment_t& env);

	void setGaugeAction(GaugeAction* _gaugeAction);
	GaugeAction* getGaugeAction() const;

	void setFermionAction(NFlavorFermionAction* _fermionAction);
	NFlavorFermionAction* getFermionAction() const;
private:
	//The gluon part of the action
	GaugeAction* gaugeAction;
	//The fermion part of the action
	NFlavorFermionAction* fermionAction;
};

} /* namespace Update */
#endif /* NFLAVORQCDACTION_H_ */
