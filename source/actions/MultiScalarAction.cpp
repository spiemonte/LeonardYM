#include "MultiScalarAction.h"

namespace Update {

MultiScalarAction::MultiScalarAction() : ScalarAction(0.,0.,0.) { }

MultiScalarAction::~MultiScalarAction() {
	for (unsigned int i = 0; i < flavorActions.size(); ++i) {
		delete flavorActions[i];
	}
}

long_real_t MultiScalarAction::energy(const environment_t& env) {
	long_real_t result = 0.;

	for (unsigned int i = 0; i < flavorActions.size(); ++i) {
		result += flavorActions[i]->energy(env);
	}

	return result;
}

GaugeGroup MultiScalarAction::force(const environment_t& env, int site, int mu) const {
	GaugeGroup result;
	set_to_zero(result);
	for (unsigned int i = 0; i < flavorActions.size(); ++i) {
		result += flavorActions[i]->force(env, site, mu);
	}
	
	return result;
}

void MultiScalarAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	if (flavorActions.size() == 1) {
		flavorActions[0]->updateForce(forceLattice, env);
	}
	else {
		extended_gauge_lattice_t tmp;
		for (unsigned int i = 0; i < flavorActions.size(); ++i) {
			//Initializations to zero of the force done inside updateForce
			flavorActions[i]->updateForce(tmp, env);
			if (i == 0) forceLattice = tmp;
			else {
#pragma omp parallel for
        			for (int site = 0; site < forceLattice.localsize; ++site) {
                			for (unsigned int mu = 0; mu < 4; ++mu) {
                        			forceLattice[site][mu] += tmp[site][mu];
					}
				}
			}
                }
        }
}

void MultiScalarAction::addAction(ScalarAction* action) {
	flavorActions.push_back(action);
}

}

