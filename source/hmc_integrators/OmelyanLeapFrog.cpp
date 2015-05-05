#include "OmelyanLeapFrog.h"

namespace Update {

OmelyanLeapFrog::OmelyanLeapFrog() { }

OmelyanLeapFrog::~OmelyanLeapFrog() { }

void OmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t t_length) {
	//The step is chosen accordingly to the number of steps
	int step = t_length/numberSteps;
	real_t lambda = 0.19318332275037863;
	
	//Calculate the force
	force->updateForce(forceLattice, env);

	//Initial update of the momenta
	this->updateMomenta(momenta, forceLattice, lambda*step);

	//Update the linkConfiguration
	this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
	env.gaugeLinkConfiguration.updateHalo();
	env.synchronize();

	//Calculate the force
	force->updateForce(forceLattice, env);

	//Second Update of the momenta
	this->updateMomenta(momenta, forceLattice, (1.-2.*lambda)*step);

	//Update the linkConfiguration
	this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
	env.gaugeLinkConfiguration.updateHalo();
	env.synchronize();

	//Central update of leapfrog
	for (unsigned int i = 0; i < numberSteps-1; ++i) {
		//Calculate the force
		force->updateForce(forceLattice, env);

		//Initial update of the momenta
		this->updateMomenta(momenta, forceLattice, 2.*lambda*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Second Update of the momenta
		this->updateMomenta(momenta, forceLattice, (1.-2.*lambda)*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();
	}

	//Calculate the force
	force->updateForce(forceLattice, env);

	//Final update of the momenta
	this->updateMomenta(momenta, forceLattice, lambda*step);
}

void OmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) {
	this->integrate(env, momenta, force, numberSteps, t_length, 0);
}

void OmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length, int forceIndex) {
	//The step is chosen accordingly to the number of steps
	real_t step = t_length/numberSteps[forceIndex];
	real_t lambda = 0.19318332275037863;
	
	//Calculate the force
	force[forceIndex]->updateForce(forceLattice, env);

	//Initial update of the momenta
	this->updateMomenta(momenta, forceLattice, lambda*step);

	if (forceIndex + 1 == force.size()) {
		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();
	} else {
		//Nested integrator, we suppose that the next force must be integrated with more precision
		this->integrate(env, momenta, force, numberSteps, step/2., forceIndex+1);
	}

	//Calculate the force
	force[forceIndex]->updateForce(forceLattice, env);

	//Second Update of the momenta
	this->updateMomenta(momenta, forceLattice, (1.-2.*lambda)*step);

	//Update the linkConfiguration
	if (forceIndex + 1 == force.size()) {
		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();
	} else {
		//Nested integrator, we suppose that the next force must be integrated with more precision
		this->integrate(env, momenta, force, numberSteps, step/2., forceIndex+1);
	}

	//Central update of leapfrog
	for (unsigned int i = 0; i < numberSteps[forceIndex]-1; ++i) {
		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Initial update of the momenta
		this->updateMomenta(momenta, forceLattice, 2.*lambda*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, step/2., forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Second Update of the momenta
		this->updateMomenta(momenta, forceLattice, (1.-2.*lambda)*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, step/2.);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, step/2., forceIndex+1);
		}
	}

	//Calculate the force
	force[forceIndex]->updateForce(forceLattice, env);

	//Final update of the momenta
	this->updateMomenta(momenta, forceLattice, lambda*step);
}

}
