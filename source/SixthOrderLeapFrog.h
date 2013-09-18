/*
 * SixthOrderLeapFrog.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef SIXTHORDERLEAPFROG_H_
#define SIXTHORDERLEAPFROG_H_

#include "Integrate.h"
#include "FourthOrderLeapFrog.h"

namespace Update {

class SixthOrderLeapFrog: public Update::Integrate {
public:
	SixthOrderLeapFrog();
	~SixthOrderLeapFrog();

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step);
	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length);
private:
	FourthOrderLeapFrog fourthOrderLeapFrog;
};

} /* namespace Update */
#endif /* LEAPFROG_H_ */
