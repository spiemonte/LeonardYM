/*
 * LeapFrog.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#include "LeapFrog.h"

namespace Update {

LeapFrog::LeapFrog() { }

LeapFrog::~LeapFrog() { }

void LeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step) {
	//Calculate the force
	force->updateForce(forceLattice, env);

	//Initial update of the momenta
	this->updateMomenta(momenta, forceLattice, step/2.);

	//Central update of leapfrog
	for (unsigned int i = 0; i < numberSteps-1; ++i) {
		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update the momenta
		this->updateMomenta(momenta, forceLattice, step);
	}

	//Final Update of the linkConfiguration
	this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step);
	env.gaugeLinkConfiguration.updateHalo();
	env.synchronize();

	//Calculate the force
	force->updateForce(forceLattice, env);

	//Final update of the momenta
	this->updateMomenta(momenta, forceLattice, step/2.);

}

void LeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) {
	this->integrate(env, momenta, force, numberSteps, t_length, 0);
}

void LeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length, int forceIndex) {
	//The step is chosen accordingly to the number of steps
	real_t step = t_length/numberSteps[forceIndex];

	//Calculate the force
	force[forceIndex]->updateForce(forceLattice, env);

	//Initial update of the momenta
	this->updateMomenta(momenta, forceLattice, step/2.);

	//Central update of leapfrog
	for (unsigned int i = 0; i < numberSteps[forceIndex]-1; ++i) {
		if (forceIndex + 1 == static_cast<int>(force.size())) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, step, forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update the momenta
		this->updateMomenta(momenta, forceLattice, step);
	}

	if (forceIndex + 1 == static_cast<int>(force.size())) {
		//Final Update of the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();
	} else {
		//Nested integrator, we suppose that the next force must be integrated with more precision
		this->integrate(env, momenta, force, numberSteps, step, forceIndex+1);
	}

	//Calculate the force
	force[forceIndex]->updateForce(forceLattice, env);

	//Final update of the momenta
	this->updateMomenta(momenta, forceLattice, step/2.);
}

} /* namespace Update */
