/*
 * PureGaugeOverrelaxation.cpp
 *
 *  Created on: May 25, 2012
 *      Author: spiem_01
 */

#include "PureGaugeOverrelaxation.h"
#include "Checkerboard.h"

namespace Update {

#if NUMCOLORS == 2
PureGaugeOverrelaxation::PureGaugeOverrelaxation() : LatticeSweep() { }
#endif
#if NUMCOLORS > 2
#ifndef MULTITHREADING
PureGaugeOverrelaxation::PureGaugeOverrelaxation() : LatticeSweep(), acceptance(0), nsteps(0),  randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
PureGaugeOverrelaxation::PureGaugeOverrelaxation() : LatticeSweep(), acceptance(0), nsteps(0) {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (unsigned int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#endif

#if NUMCOLORS == 2
PureGaugeOverrelaxation::PureGaugeOverrelaxation(const PureGaugeOverrelaxation& copy) : LatticeSweep(copy) { }
#endif
#if NUMCOLORS > 2
#ifndef MULTITHREADING
PureGaugeOverrelaxation::PureGaugeOverrelaxation(const PureGaugeOverrelaxation& copy) : LatticeSweep(), acceptance(copy.acceptance), nsteps(copy.nsteps),  randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
PureGaugeOverrelaxation::PureGaugeOverrelaxation(const PureGaugeOverrelaxation& copy) : LatticeSweep(), acceptance(copy.acceptance), nsteps(copy.nsteps) {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#endif

#if NUMCOLORS == 2
PureGaugeOverrelaxation::~PureGaugeOverrelaxation() { }
#endif
#if NUMCOLORS > 2
#ifndef MULTITHREADING
PureGaugeOverrelaxation::~PureGaugeOverrelaxation() {
	std::cout << "Overrelaxation acceptance: " << static_cast<double>(acceptance)/nsteps << std::endl;
}
#endif
#ifdef MULTITHREADING
PureGaugeOverrelaxation::~PureGaugeOverrelaxation() {
	std::cout << "Overrelaxation acceptance: " << static_cast<double>(acceptance)/nsteps << std::endl;
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomUniform[i];
	}
	delete[] randomGenerator;
	delete[] randomUniform;
}
#endif
#endif

void PureGaugeOverrelaxation::execute(environment_t& environment) {
	//Get the gauge action
	GaugeAction* gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

#ifdef MULTITHREADING
	for (unsigned int color = 0; color < checkerboard->getNumberLoops(); ++color) {
		#pragma omp parallel for //shared(environment, color) firstprivate(gaugeAction, checkerboard) default(none) schedule(dynamic)
#endif
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
#if NUMCOLORS > 2
					//take the staple
					GaugeGroup staple = gaugeAction->staple(environment.gaugeLinkConfiguration, site, mu);
					GaugeGroup Q, R;
					qr(Q,R,staple);//Compute the QR decomposition
					//Now set the det(Q) to 1:
					std::complex<real_t> determinant = conj(det(Q));
					// by multiply the first row with the complex conjugate of det:
					for (unsigned int i = 0; i < numberColors; ++i) {
						Q(0,i) *= determinant;
					}
					GaugeGroup trial = htrans(Q)*htrans(environment.gaugeLinkConfiguration[site][mu])*htrans(Q);
					real_t delta = gaugeAction->deltaAction(environment.gaugeLinkConfiguration, trial, staple, site, mu);
					//Do the accept/reject metropolis
					if (delta < 0.) {
						++acceptance;
						environment.gaugeLinkConfiguration[site][mu] = trial;
					}
#ifndef MULTITHREADING
					else if (randomUniform() < exp(-delta)) {
#endif
#ifdef MULTITHREADING
					else if ((*randomUniform[omp_get_thread_num()])() < exp(-delta)) {
#endif
						++acceptance;
						environment.gaugeLinkConfiguration[site][mu] = trial;
					}
					++nsteps;
#endif
#if NUMCOLORS == 2
					//take the staple
					GaugeGroup staple = gaugeAction->staple(environment.gaugeLinkConfiguration, site, mu);
					real_t detStaple = abs(det(staple));
					//U_new = S^\dag U^\dag S^\dag
					environment.gaugeLinkConfiguration[site][mu] = (htrans(staple)*htrans(environment.gaugeLinkConfiguration[site][mu])*htrans(staple));
					//normalize the determinant
					for (unsigned int v = 0; v < 2; ++v) {
						for (unsigned int u = 0; u < 2; ++u) {
							environment.gaugeLinkConfiguration[site][mu].at(v,u) /= detStaple;
						}
					}
#endif
#ifdef MULTITHREADING
				}
#endif
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif
	environment.synchronize();
	delete gaugeAction;

}

} /* namespace Update */
