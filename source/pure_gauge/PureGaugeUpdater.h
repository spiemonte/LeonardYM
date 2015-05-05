/*
 * PureGaugeUpdater.h
 *
 *  Created on: Feb 27, 2012
 *      Author: spiem_01
 */

#ifndef PUREGAUGEUPDATER_H_
#define PUREGAUGEUPDATER_H_

#include "LatticeSweep.h"
#include "utils/RandomSeed.h"
#include "actions/GaugeAction.h"

namespace Update {

class PureGaugeUpdater: public Update::LatticeSweep {
public:
	PureGaugeUpdater();
	PureGaugeUpdater(const PureGaugeUpdater&);
	~PureGaugeUpdater();

	/**
	 * This function perform a pure gauge update on the whole (fundamental) lattice, using the standard Wilson action
	 */
	void execute(environment_t& environment);

public:
	void updateLink(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action, double beta);

private:
#ifdef MULTITHREADING
	//The generator of random numbers
	random_generator_t** randomGenerator;
	//The generator of random numbers, uniform distribution [0,1]
	random_uniform_generator_t** randomUniform;
#endif
#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The generator of random numbers, uniform distribution [0,1]
	random_uniform_generator_t randomUniform;
#endif

	/**
	 * This function generates the radius of the SU(2) heatbath matrix using the Kennedy-Pendleton algorithm
	 * \param b the inverse of beta (Read)
	 * \return the radius r
	 */
	real_t generate_radius(real_t b);

	/**
	 * This function generates a random vector uniformly distributed on a sphere of radius radius.
	 * \param radius the radius of the sphere (Read)
	 * \param u1 the first component of the vector (Write)
	 * \param u2 the second component of the vector (Write)
	 * \param u3 the third component of the vector (Write)
	 */
	void generate_vector(real_t radius, real_t& u1, real_t& u2, real_t& u3);
};

} /* namespace Update */

#endif /* PUREGAUGEUPDATER_H_ */
