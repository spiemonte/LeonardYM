/*
 * TestForce.h
 *
 *  Created on: Jun 4, 2012
 *      Author: spiem_01
 */

#ifndef TESTFORCE_H_
#define TESTFORCE_H_
#include "Environment.h"
#include "Energy.h"

namespace Update {

class TestForce {
public:
	TestForce();
	~TestForce();

	void testForce(const environment_t& env, Energy* energy, const GaugeGroup& force, int site, int mu);

	void genericTestForce(const environment_t& env, Energy* energy, const GaugeGroup& force, int site, int mu);

private:
	environment_t tmp_env;
};

} /* namespace Update */
#endif /* TESTFORCE_H_ */
