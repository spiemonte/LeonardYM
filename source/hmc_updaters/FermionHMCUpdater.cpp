/*
 * FermionHCMUpdater.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "FermionHMCUpdater.h"
#include "utils/RandomSeed.h"
#include <omp.h>

namespace Update {

FermionHMCUpdater::FermionHMCUpdater()
#ifndef MULTITHREADING
: HMCUpdater(),  randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))) { }
#endif
#ifdef MULTITHREADING
: HMCUpdater() {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
	}
}
#endif

FermionHMCUpdater::~FermionHMCUpdater() {
#ifdef MULTITHREADING
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomNormal[i];
	}
	delete[] randomGenerator;
	delete[] randomNormal;
#endif
}

void FermionHMCUpdater::generateGaussianDiracVector(extended_dirac_vector_t& vector) {
#pragma omp parallel for
	for (int site = 0; site < vector.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int i = 0; i < diracVectorLength; ++i) {
#ifndef MULTITHREADING
				real_t realPart = randomNormal();
				real_t imagPart = randomNormal();
#endif
#ifdef MULTITHREADING
				real_t realPart = (*randomNormal[omp_get_thread_num()])();
				real_t imagPart = (*randomNormal[omp_get_thread_num()])();
#endif
				vector[site][mu](i) = std::complex<real_t>(realPart,imagPart);
			}
		}
	}
	vector.updateHalo();
}

} /* namespace Update */
