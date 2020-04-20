#ifndef LEAPFROG_H_
#define LEAPFROG_H_

#include "Integrate.h"

namespace Update {

class LeapFrog: public Update::Integrate {
public:
	LeapFrog();
	~LeapFrog();

	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step);
	virtual void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length);
private:
	extended_gauge_lattice_t forceLattice;

	void integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length, int forceIndex);
};

} /* namespace Update */
#endif /* LEAPFROG_H_ */
