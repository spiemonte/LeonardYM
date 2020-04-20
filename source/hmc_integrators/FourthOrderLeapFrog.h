#ifndef FOURTHORDERLEAPFROG_H_
#define FOURTHORDERLEAPFROG_H_

#include "Integrate.h"
#include "LeapFrog.h"

namespace Update {

class FourthOrderLeapFrog: public Update::Integrate {
public:
	FourthOrderLeapFrog();
	~FourthOrderLeapFrog();

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step);
	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length);
private:
	LeapFrog leapFrog;
};

} /* namespace Update */
#endif /* LEAPFROG_H_ */
