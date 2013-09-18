/*
 * ReUnit.h
 *
 *  Created on: Jan 10, 2013
 *      Author: spiem_01
 */

#ifndef REUNIT_H_
#define REUNIT_H_

#include "LatticeSweep.h"

namespace Update {

class ReUnit : public LatticeSweep {
public:
	ReUnit();
	~ReUnit();

	virtual void execute(environment_t& environment);

	static double testUnitarity(const extended_gauge_lattice_t& lattice);
};

} /* namespace Update */
#endif /* REUNIT_H_ */
