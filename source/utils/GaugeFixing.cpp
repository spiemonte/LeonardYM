#include "GaugeFixing.h"
#include "ExpMap.h"

namespace Update {

GaugeFixing::GaugeFixing()
#ifndef MULTITHREADING
	: randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator)), randomGeneratorUniform(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGeneratorUniform)) { }
#endif
#ifdef MULTITHREADING
	: randomGeneratorUniform(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGeneratorUniform)) {
		randomGenerator = new random_generator_t*[omp_get_max_threads()];
		randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
			randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i]));
		}
	}
#endif

GaugeFixing::GaugeFixing(const GaugeFixing&)
#ifndef MULTITHREADING
	: randomGenerator(RandomSeed::randomSeed()), randomNormal(RandomSeed::getNormalNumberGenerator(randomGenerator)), randomGeneratorUniform(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGeneratorUniform)) { }
#endif
#ifdef MULTITHREADING
	: randomGeneratorUniform(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGeneratorUniform)) {
		randomGenerator = new random_generator_t*[omp_get_max_threads()];
		randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
			randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i]));
		}
	}
#endif

GaugeFixing::~GaugeFixing() {
#ifdef MULTITHREADING
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomNormal[i];
	}
	delete[] randomGenerator;
	delete[] randomNormal;
#endif
}

void GaugeFixing::generateRandomGaugeTransformation(extended_matrix_lattice_t& gauge_transformation, real_t epsilon) {
	 ExponentialMap expMap;

#pragma omp parallel for
	for (int site = 0; site < gauge_transformation.localsize; ++site) {
		//Antihermitian part
		for (int i = 0; i < numberColors; ++i) {
			for (int j = i+1; j < numberColors; ++j) {
#ifndef MULTITHREADING
				real_t realPart = epsilon*randomNormal()/2.;
				real_t imagPart = epsilon*randomNormal()/2.;
#endif
#ifdef MULTITHREADING
				real_t realPart = epsilon*(*randomNormal[omp_get_thread_num()])()/2.;
				real_t imagPart = epsilon*(*randomNormal[omp_get_thread_num()])()/2.;
#endif
				gauge_transformation[site].at(i,j) = std::complex<real_t>(realPart,imagPart);
				gauge_transformation[site].at(j,i) = std::complex<real_t>(-realPart,imagPart);
			}
		}
		for (int i = 0; i < numberColors; ++i) {
			gauge_transformation[site].at(i,i) = 0.;
		}
		//Antihermitian Traceless part
		for (int i = 1; i < numberColors; ++i) {
#ifndef MULTITHREADING
			real_t imagPart = epsilon*randomNormal()/2.;
#endif
#ifdef MULTITHREADING
			real_t imagPart = (*randomNormal[omp_get_thread_num()])()/2.;
#endif
			for (int j = 0; j < i; ++j) {
				gauge_transformation[site].at(j,j) += std::complex<real_t>( 0, imagPart/sqrt(static_cast<real_t>(i*(i+1)/2.)) );
			}
			gauge_transformation[site].at(i,i) += std::complex<real_t>( 0,-imagPart*i*sqrt( static_cast<real_t>( 2./(i*(i+1)) ) ) );
		}
		
		gauge_transformation[site] = expMap.exp(gauge_transformation[site]);
	}

	gauge_transformation.updateHalo();
}

void GaugeFixing::transform(extended_gauge_lattice_t& lattice, const extended_matrix_lattice_t& gauge_transformation) const {
	typedef extended_matrix_lattice_t LT;

#pragma omp parallel for
	for (int site = 0; site < gauge_transformation.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) lattice[site][mu] = gauge_transformation[site]*lattice[site][mu]*htrans(gauge_transformation[LT::sup(site,mu)]);
	}

	lattice.updateHalo();
}

void GaugeFixing::getLieAlgebraField(extended_gauge_lattice_t& Afield, const extended_gauge_lattice_t& lattice) const {
	typedef extended_matrix_lattice_t LT;

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			Afield[site][mu] = lattice[site][mu] - htrans(lattice[site][mu]);
			Afield[site][mu] = Afield[site][mu]/(std::complex<real_t>(0.,2.));
			std::complex<real_t> tr = trace(Afield[site][mu]);
			for (int i = 0; i < numberColors; ++i) {
				Afield[site][mu].at(i,i) -= tr/static_cast<real_t>(numberColors);
			}
		}
	}

	Afield.updateHalo();
}

bool GaugeFixing::metropolis(real_t delta) {
        if (delta < -50.) return false;
        else if (delta >= 0.) return true;
	
	real_t rU = 0;
	if (isOutputProcess()) {
		rU = randomUniform();
	}
	reduceAllSum(rU);

        if (exp(delta) > rU) return true;
        else return false;
}

}
