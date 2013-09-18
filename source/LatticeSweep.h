/*
 * LatticeSweep.h
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#ifndef LATTICESWEEP_H_
#define LATTICESWEEP_H_
#include "Environment.h"

namespace Update {

class LatticeSweep {
public:
	LatticeSweep();
	LatticeSweep(unsigned int _numberTimes, unsigned int _sweepToJump);
	virtual ~LatticeSweep();

	static LatticeSweep* getInstance(const std::string& name);
	static LatticeSweep* read(const std::string& toRead);

	void call(environment_t& environment);

	virtual void execute(environment_t& environment) = 0;

	void setNumberTimes(unsigned int _numberTimes);
	void setSweepToJump(unsigned int _sweepToJump);

	static void printSweepsName();
private:
	unsigned int numberTimes;
	unsigned int sweepToJump;
};

} /* namespace Update */
#endif /* LATTICESWEEP_H_ */
