/*
 * WilsonLoop.h
 *
 *  Created on: Jul 30, 2012
 *      Author: spiem_01
 */

#ifndef WILSONLOOP_H_
#define WILSONLOOP_H_
#include "LatticeSweep.h"

namespace Update {

class WilsonLoop : public LatticeSweep {
public:
	WilsonLoop();
	~WilsonLoop();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* WILSONLOOP_H_ */
