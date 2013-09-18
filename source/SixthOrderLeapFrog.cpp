#include "SixthOrderLeapFrog.h"

namespace Update {

SixthOrderLeapFrog::SixthOrderLeapFrog() { }

SixthOrderLeapFrog::~SixthOrderLeapFrog() { }

void SixthOrderLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step) {
	real_t b1 = 1./(2.-pow(2.,1./5.));
	real_t b2 = -pow(2.,1./5.)/(2.-pow(2.,1./5.));
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b1*step);
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b2*step);
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b1*step);
}

void SixthOrderLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) {
	real_t b1 = 1./(2.-pow(2.,1./5.));
	real_t b2 = -pow(2.,1./5.)/(2.-pow(2.,1./5.));
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b1*t_length);
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b2*t_length);
	fourthOrderLeapFrog.integrate(env, momenta, force, numberSteps, b1*t_length);
}

}
