#ifndef MULTISCALARACTION_H
#define MULTISCALARACTION_H
#include "Environment.h"
#include "ScalarAction.h"

namespace Update {

class MultiScalarAction : public ScalarAction {
public:
	MultiScalarAction();
	~MultiScalarAction();

#ifndef __IBMCPP__
        using Force::force;
#endif
	
	void addAction(ScalarAction*);
	
	virtual long_real_t energy(const environment_t& env);
	virtual GaugeGroup force(const environment_t& env, int site, int mu) const;
    virtual void updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env);
private:
    std::vector<ScalarAction*> flavorActions;

	MultiScalarAction(const MultiScalarAction&) : ScalarAction() { }
};

} /* namespace Update */

#endif /* ADJOINTSCALARACTION_H_ */

