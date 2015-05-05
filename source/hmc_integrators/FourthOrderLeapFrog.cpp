#include "FourthOrderLeapFrog.h"

namespace Update {

FourthOrderLeapFrog::FourthOrderLeapFrog() { }

FourthOrderLeapFrog::~FourthOrderLeapFrog() { }

void FourthOrderLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, Force* force, int numberSteps, real_t step) {
	real_t b1 = 1./(2.-pow(2.,1./3.));
	real_t b2 = -pow(2.,1./3.)/(2.-pow(2.,1./3.));
	leapFrog.integrate(env, momenta, force, numberSteps, b1*step);
	leapFrog.integrate(env, momenta, force, numberSteps, b2*step);
	leapFrog.integrate(env, momenta, force, numberSteps, b1*step);
}

void FourthOrderLeapFrog::integrate(environment_t& env, extended_gauge_lattice_t& momenta, const std::vector<Force*>& force, const std::vector<unsigned int>& numberSteps, real_t t_length) {
	real_t b1 = 1./(2.-pow(2.,1./3.));
	real_t b2 = -pow(2.,1./3.)/(2.-pow(2.,1./3.));
	leapFrog.integrate(env, momenta, force, numberSteps, b1*t_length);
	leapFrog.integrate(env, momenta, force, numberSteps, b2*t_length);
	leapFrog.integrate(env, momenta, force, numberSteps, b1*t_length);
}


}
