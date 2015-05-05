/*
 * Plaquette.h
 *
 *  Created on: Feb 29, 2012
 *      Author: spiem_01
 */

#ifndef PLAQUETTE_H_
#define PLAQUETTE_H_

#include "LatticeSweep.h"

namespace Update {

class Plaquette: public Update::LatticeSweep {
public:
	Plaquette();
	~Plaquette();

	/**
	 * This function measure the plaquette on the lattice and it prints out its value.
	 * @param enviroment
	 * @param sweep
	 * @param n
	 */
	virtual void execute(environment_t& environment);

	static long_real_t temporalPlaquette(const extended_gauge_lattice_t& gaugeLinkConfiguration);
};

} /* namespace Update */
#endif /* PLAQUETTE_H_ */
