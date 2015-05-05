/*
 * StochasticEstimator.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

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
	void generateRandomNoise(extended_dirac_vector_t& vector);

	void generateRandomNoise(extended_dirac_vector_t* vector, int t);

	void generateSource(extended_dirac_vector_t& vector, int alpha, int c);

	template<typename T> T mean(const std::vector<T>& v) {
		typename std::vector<T>::const_iterator i;
		T res(0.);
		for (i = v.begin(); i != v.end(); ++i) {
			res += *i;
		}
		return res/v.size();
	}

	template<typename T> T standardDeviation(const std::vector<T>& v) {
		T mn(mean(v));
		T res(0.);
		typename std::vector<T>::const_iterator i;
		for (i = v.begin(); i != v.end(); ++i) {
			res += ((*i) - mn)*((*i) - mn);
		}
		return sqrt(res/(v.size()-1));
	}

private:
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
