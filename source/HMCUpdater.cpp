/*
 * HMCUpdater.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#include "HMCUpdater.h"
#include <omp.h>

namespace Update {

HMCUpdater::HMCUpdater()
#ifndef MULTITHREADING
	: randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator)), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)), acceptance(0), counter(0) { }
#endif
#ifdef MULTITHREADING
	: singleRandomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(singleRandomGenerator)), acceptance(0), counter(0) {
		randomGenerator = new random_generator_t*[omp_get_max_threads()];
		randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
			randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i]));
		}
	}
#endif

HMCUpdater::~HMCUpdater() {
	if (isOutputProcess()) std::cout << "Acceptance rate: " << static_cast<double>(acceptance)/counter << std::endl;
#ifdef MULTITHREADING
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomNormal[i];
	}
	delete[] randomGenerator;
	delete[] randomNormal;
#endif
}

void HMCUpdater::randomMomenta(extended_gauge_lattice_t& momenta) {
#pragma omp parallel for
	for (int position = 0; position < momenta.localsize; ++position) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			//Antihermitian part
			for (int i = 0; i < numberColors; ++i) {
				for (int j = i+1; j < numberColors; ++j) {
#ifndef MULTITHREADING
					real_t realPart = randomNormal()/2.;
					real_t imagPart = randomNormal()/2.;
#endif
#ifdef MULTITHREADING
					real_t realPart = (*randomNormal[omp_get_thread_num()])()/2.;
					real_t imagPart = (*randomNormal[omp_get_thread_num()])()/2.;
#endif
					momenta[position][mu].at(i,j) = std::complex<real_t>(realPart,imagPart);
					momenta[position][mu].at(j,i) = std::complex<real_t>(-realPart,imagPart);
				}
			}
			for (int i = 0; i < numberColors; ++i) {
				momenta[position][mu].at(i,i) = 0.;
			}
			//Antihermitian Traceless part
			for (int i = 1; i < numberColors; ++i) {
#ifndef MULTITHREADING
				real_t imagPart = randomNormal()/2.;
#endif
#ifdef MULTITHREADING
				real_t imagPart = (*randomNormal[omp_get_thread_num()])()/2.;
#endif
				for (int j = 0; j < i; ++j) {
					momenta[position][mu].at(j,j) += std::complex<real_t>( 0, imagPart/sqrt(static_cast<real_t>(i*(i+1)/2.)) );
				}
				momenta[position][mu].at(i,i) += std::complex<real_t>( 0,-imagPart*i*sqrt( static_cast<real_t>( 2./(i*(i+1)) ) ) );
			}
		}
	}
	momenta.updateHalo();//TODO is needed?
}

bool HMCUpdater::metropolis(long_real_t energyOld, long_real_t energyNew) {
	int decision = 0;
	if (isOutputProcess()) {
		++counter;
		long_real_t delta = energyNew - energyOld;
		std::cout << "Delta energy in metropolis: " << delta;
		if (delta < 0.) {
			++acceptance;
			decision = 1;
		}
		else if (randomUniform() < exp(-delta)) {
			++acceptance;
			decision = 1;
		}
		else decision = 0;
	}
#ifdef ENABLE_MPI
	MPI_Bcast(&decision, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
	if (decision == 0) {
		if (isOutputProcess()) std::cout << " - rejected!" << std::endl;
		return false;
	}
	else {
		if (isOutputProcess()) std::cout << " - accepted!" << std::endl;
		return true;
	}
}

bool HMCUpdater::metropolis(long_real_t value) {
	int decision = 0;
        if (isOutputProcess()) {
                ++counter;
                std::cout << "Metropolis distribution: min(1," << value << ") ";
                if (value > 1.) {
                        ++acceptance;
                        decision = 1;
                }
                else if (randomUniform() < value) {
                        ++acceptance;
                        decision = 1;
                }
                else decision = 0;
        }
#ifdef ENABLE_MPI
        MPI_Bcast(&decision, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
        if (decision == 0) {
                if (isOutputProcess()) std::cout << " - rejected!" << std::endl;
                return false;
        }
        else {
                if (isOutputProcess()) std::cout << " - accepted!" << std::endl;
                return true;
        }
}

long_real_t HMCUpdater::momentaEnergy(const extended_gauge_lattice_t& momenta) const {
	long_real_t energy = 0.;
#pragma omp parallel for reduction(+:energy)
	for (int position = 0; position < momenta.localsize; ++position) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			energy += -real(trace(momenta[position][mu]*momenta[position][mu]));
		}
	}
	reduceAllSum(energy);
	return energy;
}

} /* namespace Update */
