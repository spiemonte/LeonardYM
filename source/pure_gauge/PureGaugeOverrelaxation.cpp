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
	typedef extended_gauge_lattice_t::Layout Layout;
	//Get the gauge action
	GaugeAction* gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

#ifndef ENABLE_MPI
#ifdef MULTITHREADING
	for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#endif
#pragma omp parallel for //shared(environment, color) firstprivate(gaugeAction, checkerboard) default(none) schedule(dynamic)
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
					this->updateLink(environment.gaugeLinkConfiguration,site,mu,gaugeAction);
#ifdef MULTITHREADING
				}
#endif
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif
#endif

	//We suppose that always MPI+MTH
#ifdef ENABLE_MPI
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
			if (processor == Layout::this_processor) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
				for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, gaugeAction);
						}
					}
				}
			}
			else {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
				for (int site = environment.gaugeLinkConfiguration.sharedsize; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, gaugeAction);
						}
					}
				}
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
	}
#endif
	
	environment.synchronize();
	delete gaugeAction;

}

void PureGaugeOverrelaxation::updateLink(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* gaugeAction) {
#if NUMCOLORS > 2
	//take the staple
	GaugeGroup staple = gaugeAction->staple(lattice, site, mu);
	GaugeGroup Q, R;
	qr(Q,R,staple);//Compute the QR decomposition
	//Now set the det(Q) to 1:
	std::complex<real_t> determinant = conj(det(Q));
	// by multiply the first row with the complex conjugate of det:
	for (unsigned int i = 0; i < numberColors; ++i) {
		Q(0,i) *= determinant;
	}
	GaugeGroup trial = htrans(Q)*htrans(lattice[site][mu])*htrans(Q);
	real_t delta = gaugeAction->deltaAction(lattice, trial, staple, site, mu);
	//Do the accept/reject metropolis
	if (delta < 0.) {
		++acceptance;
		lattice[site][mu] = trial;
	}
#ifndef MULTITHREADING
	else if (randomUniform() < exp(-delta)) {
#endif
#ifdef MULTITHREADING
	else if ((*randomUniform[omp_get_thread_num()])() < exp(-delta)) {
#endif
		++acceptance;
		lattice[site][mu] = trial;
	}
	++nsteps;
#endif
#if NUMCOLORS == 2
	//take the staple
	GaugeGroup staple = gaugeAction->staple(lattice, site, mu);
	real_t detStaple = abs(det(staple));
	//U_new = S^\dag U^\dag S^\dag
	lattice[site][mu] = (htrans(staple)*htrans(lattice[site][mu])*htrans(staple));
	//normalize the determinant
	for (unsigned int v = 0; v < 2; ++v) {
		for (unsigned int u = 0; u < 2; ++u) {
			lattice[site][mu].at(v,u) /= detStaple;
		}
	}
#endif
}

} /* namespace Update */
