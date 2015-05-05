/*
 * Integrate.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include "hmc_forces/Force.h"
#include <string>

namespace Update {

class Integrate {
public:
	Integrate();
	virtual ~Integrate();

	static Integrate* getInstance(const std::string& nameAlgorithm);

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step) = 0;

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) = 0;

	void updateLinkConfiguration(extended_gauge_lattice_t& linkConfiguration, const extended_gauge_lattice_t& momenta, real_t epsilon);
	void updateMomenta(extended_gauge_lattice_t& momenta, const extended_gauge_lattice_t& force, real_t epsilon);
};

} /* namespace Update */
#endif /* INTEGRATE_H_ */
