#include "RandomScalarUpdater.h"
#include "utils/RandomSeed.h"
#include <omp.h>

namespace Update {

RandomScalarUpdater::RandomScalarUpdater()
#ifndef MULTITHREADING
: randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))) { }
#endif
#ifdef MULTITHREADING
{
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
	}
}
#endif

RandomScalarUpdater::RandomScalarUpdater(const RandomScalarUpdater&)
#ifndef MULTITHREADING
: LatticeSweep(), randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator,sqrt(0.5))) { }
#endif
#ifdef MULTITHREADING
: LatticeSweep()
{
        randomGenerator = new random_generator_t*[omp_get_max_threads()];
        randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
        for (int i = 0; i < omp_get_max_threads(); ++i) {
                randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
                randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],sqrt(0.5)));
        }
}
#endif

RandomScalarUpdater::~RandomScalarUpdater() {
#ifdef MULTITHREADING
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomNormal[i];
	}
	delete[] randomGenerator;
	delete[] randomNormal;
#endif
}

void RandomScalarUpdater::execute(environment_t& environment) {
	unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");
        if (environment.adjoint_scalar_fields.size() != aNf) {
                if (isOutputProcess()) std::cout << "RandomScalarUpdater::Setting the adjoint scalar Nf to " << aNf << std::endl;
                environment.adjoint_scalar_fields.resize(aNf);
	}
        
	for (unsigned int f = 0; f < aNf; ++f) {
                this->generateGaussianAdjointScalar(environment.adjoint_scalar_fields[f]);
        }

	unsigned int nf = environment.configurations.get<unsigned int>("fundamental_nf_scalars");
        if (environment.fundamental_scalar_fields.size() != nf) {
                if (isOutputProcess()) std::cout << "RandomScalarUpdater::Setting the fundamental scalar Nf to " << nf << std::endl;
                environment.fundamental_scalar_fields.resize(nf);
        }

        for (unsigned int f = 0; f < nf; ++f) {
                this->generateGaussianFundamentalScalar(environment.fundamental_scalar_fields[f]);
        }
}

void RandomScalarUpdater::generateGaussianAdjointScalar(extended_adjoint_real_color_vector_t& vector) {
#pragma omp parallel for
	for (int site = 0; site < vector.localsize; ++site) {
		for (unsigned int i = 0; i < numberColors*numberColors - 1; ++i) {
#ifndef MULTITHREADING
			real_t realPart = randomNormal();
#endif
#ifdef MULTITHREADING
			real_t realPart = (*randomNormal[omp_get_thread_num()])();
#endif
			vector[site](i) = realPart;
		}
	}
	vector.updateHalo();
}

void RandomScalarUpdater::generateGaussianFundamentalScalar(extended_color_vector_t& vector) {
#pragma omp parallel for
        for (int site = 0; site < vector.localsize; ++site) {
                for (unsigned int i = 0; i < numberColors; ++i) {
#ifndef MULTITHREADING
                        real_t realPart = randomNormal();
                        real_t imagPart = randomNormal();
#endif
#ifdef MULTITHREADING
                        real_t realPart = (*randomNormal[omp_get_thread_num()])();
                        real_t imagPart = (*randomNormal[omp_get_thread_num()])();
#endif
                        vector[site](i) = std::complex<real_t>(realPart,imagPart);
                }
        }
        vector.updateHalo();
}

} /* namespace Update */
