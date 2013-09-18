/*
 * HMCUpdater.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef HMCUPDATER_H_
#define HMCUPDATER_H_
#include "Environment.h"
#include "RandomSeed.h"

namespace Update {

/**
 * This class is a container of some algorithms used in all the HMC routines
 */
class HMCUpdater {
public:
	HMCUpdater();
	~HMCUpdater();

	/**
	 * This function random initialize the momenta, usually called in the beginning of a HMC integration
	 * @param momenta
	 */
	void randomMomenta(extended_gauge_lattice_t& momenta);

	/**
	 * This function returns the energy of the momenta (ie: -real(trace(P^2)) )
	 * @param momenta
	 * @return the energy
	 */
	long_real_t momentaEnergy(const extended_gauge_lattice_t& momenta) const;

	/**
	 * This function computes the metropolis step
	 * @param energyOld
	 * @param energyNew
	 * @return true: accepted, false rejected
	 */
	bool metropolis(long_real_t energyOld, long_real_t energyNew);
	
	bool metropolis(long_real_t value);
private:
#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The generator of normal random numbers
	random_normal_generator_t randomNormal;
#endif
#ifdef MULTITHREADING
	//The generators of random numbers
	random_generator_t** randomGenerator;
	//The generators of normal random numbers
	random_normal_generator_t** randomNormal;
	//The generator of random numbers
	random_generator_t singleRandomGenerator;
#endif
	//The generator of uniform random numbers
	random_uniform_generator_t randomUniform;
	//Acceptance for the metropolis steps
	unsigned int acceptance;
	//Global counter for the metropolis steps
	unsigned int counter;
};

} /* namespace Update */
#endif /* HMCUPDATER_H_ */
