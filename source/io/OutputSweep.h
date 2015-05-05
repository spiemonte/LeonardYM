/*
 * OutputSweep.h
 *
 *  Created on: May 25, 2012
 *      Author: spiem_01
 */

#ifndef OUTPUTSWEEP_H_
#define OUTPUTSWEEP_H_

#include "LatticeSweep.h"

namespace Update {

class OutputSweep: public Update::LatticeSweep {
public:
	OutputSweep();
	~OutputSweep();

	/**
	 * This function writes down the configuration
	 */
	void execute(environment_t& environment);
private:
	static int get_latticenumber() {
		return (latticenumber);
	}

	static const int latticenumber = 0;
};

} /* namespace Update */
#endif /* OUTPUTSWEEP_H_ */
