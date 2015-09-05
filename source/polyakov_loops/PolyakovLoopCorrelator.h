/*
 * PolyakovLoop.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef POLYAKOVLOOPCORRELATOR_H_
#define POLYAKOVLOOPCORRELATOR_H_

#include "LatticeSweep.h"

namespace Update {

class PolyakovLoopCorrelator : public Update::LatticeSweep {
public:
	PolyakovLoopCorrelator();
	~PolyakovLoopCorrelator();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* POLYAKOVLOOPCORRELATOR_H_ */
