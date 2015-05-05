/*
 * TestCommunication.h
 *
 *  Created on: Feb 25, 2013
 *      Author: spiem_01
 */

#ifndef TESTCOMMUNICATION_H_
#define TESTCOMMUNICATION_H_
#include "LatticeSweep.h"

namespace Update {

class TestCommunication : public LatticeSweep {
public:
	TestCommunication();
	~TestCommunication();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* TESTCOMMUNICATION_H_ */
