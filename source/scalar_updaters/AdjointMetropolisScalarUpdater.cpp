/*
 * AdjointMetropolisScalarUpdater.cpp
 *
 *  Created on: Apr 2, 2012
 *      Author: spiem_01
 */

#include "AdjointMetropolisScalarUpdater.h"
#include "utils/RandomSeed.h"
#include "actions/AdjointScalarAction.h"
#include <omp.h>

namespace Update {

AdjointMetropolisScalarUpdater::AdjointMetropolisScalarUpdater()
#ifndef MULTITHREADING
: LatticeSweep(), randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
: LatticeSweep()
{
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif

AdjointMetropolisScalarUpdater::AdjointMetropolisScalarUpdater(const AdjointMetropolisScalarUpdater& toCopy)
#ifndef MULTITHREADING
: LatticeSweep(toCopy), randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif
#ifdef MULTITHREADING
: LatticeSweep(toCopy) 
{
        randomGenerator = new random_generator_t*[omp_get_max_threads()];
        randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
        randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
        for (int i = 0; i < omp_get_max_threads(); ++i) {
                randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
                randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
                randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
        }
}
#endif

AdjointMetropolisScalarUpdater::~AdjointMetropolisScalarUpdater() {
#ifdef MULTITHREADING
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomNormal[i];
		delete randomUniform[i];
	}
	delete[] randomGenerator;
	delete[] randomNormal;
	delete[] randomUniform;
#endif
}

void AdjointMetropolisScalarUpdater::execute(environment_t& environment) {
	unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");
	if (aNf != environment.adjoint_scalar_fields.size()) {
		if (isOutputProcess()) std::cout << "AdjointMetropolisScalarUpdater::Number of adjoint scalar incorrectly set!" << std::endl;
		exit(109);
	} 
	
	double mu = environment.configurations.get<double>("AdjointMetropolisScalarUpdater::scalar_mass");
        double lambda = environment.configurations.get<double>("AdjointMetropolisScalarUpdater::quartic_coupling");
	double lambda_8 = environment.configurations.get<double>("AdjointMetropolisScalarUpdater::adjoint_quartic_coupling");
	double epsilon = environment.configurations.get<double>("AdjointMetropolisScalarUpdater::epsilon");
        AdjointScalarAction* action = new AdjointScalarAction(mu, lambda, lambda_8);

	unsigned int trials = environment.configurations.get<unsigned int>("AdjointMetropolisScalarUpdater::number_of_hits");

	typedef extended_adjoint_color_vector_t::Layout Layout;

	int acceptance = 0;

	std::vector<extended_adjoint_color_vector_t>::iterator scalar_field;

	/*{
		AdjointVector oldphi = environment.adjoint_scalar_fields[0][5];
		AdjointVector newphi = environment.adjoint_scalar_fields[0][5];
		newphi[0] += 0.22;
		newphi[1] -= 0.13;
		AdjointVector kinetic_coupling = action->getKineticCoupling(environment.getAdjointLattice(), environment.adjoint_scalar_fields[0], 5);
		real_t deltaMet = action->deltaEnergy(kinetic_coupling, oldphi, newphi);
		real_t deltaE = -action->energy(environment);
		environment.adjoint_scalar_fields[0][5] = newphi;
		deltaE += action->energy(environment);
		if (isOutputProcess()) std::cout << "Consistency test of the action " << deltaMet << " " << deltaE << std::endl; 
	}*/

        for (scalar_field = environment.adjoint_scalar_fields.begin(); scalar_field < environment.adjoint_scalar_fields.end(); ++scalar_field) { 
		for (int block = 0; block < 2; ++block) {
#pragma omp parallel for reduction(+:acceptance)
			for (int site = 0; site < scalar_field->localsize; ++site) {
				//White/Black partitioning
				if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == block) {
					AdjointVector kinetic_coupling = action->getKineticCoupling(environment.getAdjointLattice(), *scalar_field, site);
					for (unsigned int trial = 0; trial < trials; ++trial) {
						AdjointVector proposal = (*scalar_field)[site];
						for (unsigned int c = 0; c < numberColors*numberColors - 1; ++c) {
#ifndef MULTITHREADING
							proposal[c] += epsilon*std::complex<real_t>(randomNormal(), randomNormal());
#endif
#ifdef MULTITHREADING
							proposal[c] += epsilon*std::complex<real_t>((*randomNormal[omp_get_thread_num()])(), (*randomNormal[omp_get_thread_num()])());
#endif
						}
						real_t delta = action->deltaEnergy(kinetic_coupling, (*scalar_field)[site], proposal);
						//Do the accept/reject metropolis
        					if (delta < 0.) {
					                ++acceptance;
					                (*scalar_field)[site] = proposal;
       						 }
#ifndef MULTITHREADING
        					else if (randomUniform() < exp(-delta)) {
#endif
#ifdef MULTITHREADING
        					else if ((*randomUniform[omp_get_thread_num()])() < exp(-delta)) {
#endif
                					++acceptance;
                					(*scalar_field)[site] = proposal;
        					}
					}
					
				}	
			}
			scalar_field->updateHalo();
		}
	}
	
	reduceAllSum(acceptance);

	if (isOutputProcess()) std::cout << "AdjointMetropolisScalarUpdater::Acceptance rate: " << acceptance/static_cast<double>(Layout::globalVolume*trials*aNf) << std::endl;

	delete action;
}

void AdjointMetropolisScalarUpdater::registerParameters(po::options_description& desc) {
        static bool single = true;
        if (single) desc.add_options()
                ("AdjointMetropolisScalarUpdater::scalar_mass", po::value<double>()->default_value(0.0), "set the value of the adjoint scalar mass m*m")
                ("AdjointMetropolisScalarUpdater::quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the quartic coupling lambda")
		("AdjointMetropolisScalarUpdater::adjoint_quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the adjoint quartic coupling lambda")
		("AdjointMetropolisScalarUpdater::epsilon", po::value<double>()->default_value(0.1), "set the value of epsilon for the random update")
                ("AdjointMetropolisScalarUpdater::number_of_hits", po::value<unsigned int>()->default_value(1), "set the number of the trials for the multihit metropolis")
                ;
        single = false;
}

} /* namespace Update */
