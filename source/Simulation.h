/*
 * Simulation.h
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_
#include "Environment.h"
#include "LatticeSweep.h"
#include <list>

namespace Update {

class Simulation {
public:
	Simulation(const environment_t& _environment);
	~Simulation();

	/**
	 * This function starts the lattice, the environment and all other factors for the simulation
	 */
	void starter();

	void warmUp();

	void measurement();

private:
	std::list<LatticeSweep*> listWarmUpSweeps;
	std::list<LatticeSweep*> listMeasurementSweeps;
	environment_t environment;
};

} /* namespace Update */
#endif /* SIMULATION_H_ */
