#ifndef STOCHASTICESTIMATOR_H_
#define STOCHASTICESTIMATOR_H_
#include "Environment.h"
#include "utils/RandomSeed.h"


namespace Update {

class StochasticEstimator {
public:
	StochasticEstimator();
	StochasticEstimator(const StochasticEstimator& copy);
	~StochasticEstimator();

protected:
	template<typename TVector> void generateRandomNoise(TVector& vector) {
#pragma omp parallel for
		for (int site = 0; site < vector.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int i = 0; i < diracVectorLength; ++i) {
#ifndef MULTITHREADING
					real_t realPart = (randomInteger() == 0 ? -1 : 1);
#endif
#ifdef MULTITHREADING
					real_t realPart = ((*randomInteger[omp_get_thread_num()])() == 0 ? -1 : 1);
#endif
					vector[site][mu][i] = std::complex<real_t>(realPart,0.);
				}
			}
		}
		vector.updateHalo();
	}

	void generateRandomNoise(extended_dirac_vector_t* vector, int t);

	void generateRandomNoise(extended_dirac_vector_t& vector, int t0, int t1);

	void generateSource(extended_dirac_vector_t& vector, int alpha, int c);

	void generateMomentumSource(extended_dirac_vector_t& vector, std::vector<real_t> p, int alpha, int c);

	void smearSource(extended_dirac_vector_t& vector, const extended_fermion_lattice_t& lattice, unsigned int levels, const real_t& alpha, int no_smear_dir = 3, const real_t& K = 0.9);

	template<typename T> T mean(const std::vector<T>& v) {
		typename std::vector<T>::const_iterator i;
		T res(0.);
		for (i = v.begin(); i != v.end(); ++i) {
			res += *i;
		}
		return res/static_cast<T>(v.size());
	}

	template<typename T> T standardDeviation(const std::vector<T>& v) {
		T mn(mean(v));
		T res(0.);
		typename std::vector<T>::const_iterator i;
		for (i = v.begin(); i != v.end(); ++i) {
			res += ((*i) - mn)*((*i) - mn);
		}
		return sqrt(res/static_cast<T>(v.size()-1));
	}

#ifndef MULTITHREADING
	//The generator of random numbers
	random_generator_t randomGenerator;
	//The random integer generator
	random_integer_generator_t randomInteger;
#endif
#ifdef MULTITHREADING
	//The generator of random numbers
	random_generator_t** randomGenerator;
	//The random integer generator
	random_integer_generator_t** randomInteger;
#endif
};

} /* namespace Update */
#endif /* STOCHASTICESTIMATOR_H_ */
