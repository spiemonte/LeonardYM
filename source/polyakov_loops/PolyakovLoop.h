/*
 * PolyakovLoop.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef POLYAKOVLOOP_H_
#define POLYAKOVLOOP_H_

#include "LatticeSweep.h"

namespace Update {

class PolyakovLoop: public Update::LatticeSweep {
public:
	PolyakovLoop();
	~PolyakovLoop();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* POLYAKOVLOOP_H_ */
