/*
 * FourthOmelyanLeapFrog.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef FOURTHOMELYANLEAPFROG_H_
#define FOURTHOMELYANLEAPFROG_H_

#include "Integrate.h"

namespace Update {

class FourthOmelyanLeapFrog : public Update::Integrate {
public:
	FourthOmelyanLeapFrog();
	~FourthOmelyanLeapFrog();

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step);
	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length);

private:
	extended_gauge_lattice_t forceLattice;
	
	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length, int forceIndex);
};

} /* namespace Update */
#endif /* LEAPFROG_H_ */
