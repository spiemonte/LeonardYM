#include "FourthOmelyanLeapFrog.h"

namespace Update {

FourthOmelyanLeapFrog::FourthOmelyanLeapFrog() { }

FourthOmelyanLeapFrog::~FourthOmelyanLeapFrog() { }

void FourthOmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t t_length) {
	//The step is chosen accordingly to the number of steps
	int step = t_length/numberSteps;
	real_t theta = 0.083983152628;
	real_t rho = 0.2539785108410595;
	real_t mu = -0.03230286765269967;
	real_t lambda = 0.6822365335719091;
	
	for (unsigned int i = 0; i < numberSteps; ++i) {
		//Position version of the integrator, first we integrate the linkconf
		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, theta*step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, rho*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, lambda*step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, mu*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, (1. - 2.*(lambda + theta))*step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, (1. - 2.*(mu + rho))*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, (1. - 2.*(lambda + theta))*step/2.);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, mu*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, lambda*step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();

		//Calculate the force
		force->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, rho*step);

		//Update the linkConfiguration
		this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, theta*step);
		env.gaugeLinkConfiguration.updateHalo();
		env.synchronize();
	}
}

void FourthOmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) {
	this->integrate(env, momenta, force, numberSteps, t_length, 0);
}

void FourthOmelyanLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length, int forceIndex) {
	//The step is chosen accordingly to the number of steps
	real_t step = t_length/numberSteps[forceIndex];
	real_t theta = 0.083983152628;
	real_t rho = 0.2539785108410595;
	real_t mu = -0.03230286765269967;
	real_t lambda = 0.6822365335719091;

	for (unsigned int i = 0; i < numberSteps[forceIndex]; ++i) {
		//Position version of the integrator, first we integrate the linkconf
		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, theta*step);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, theta*step, forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, rho*step);

		//Update the linkConfiguration
		if (forceIndex + 1 == force.size()) {
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, lambda*step);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, lambda*step, forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, mu*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, (1. - 2.*(lambda + theta))*step/2.);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, (1. - 2.*(lambda + theta))*step/2., forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, (1. - 2.*(mu + rho))*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, (1. - 2.*(lambda + theta))*step/2.);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, (1. - 2.*(lambda + theta))*step/2., forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, mu*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, lambda*step);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, lambda*step, forceIndex+1);
		}

		//Calculate the force
		force[forceIndex]->updateForce(forceLattice, env);

		//Update of the momenta
		this->updateMomenta(momenta, forceLattice, rho*step);

		if (forceIndex + 1 == force.size()) {
			//Update the linkConfiguration
			this->updateLinkConfiguration(env.gaugeLinkConfiguration, momenta, theta*step);
			env.gaugeLinkConfiguration.updateHalo();
			env.synchronize();
		} else {
			//Nested integrator, we suppose that the next force must be integrated with more precision
			this->integrate(env, momenta, force, numberSteps, theta*step, forceIndex+1);
		}
	}
}

}
