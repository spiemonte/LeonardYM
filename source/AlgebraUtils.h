/*
 * AlgebraUtils.h
 *
 *  Created on: Apr 20, 2012
 *      Author: spiem_01
 */

#ifndef ALGEBRAUTILS_H_
#define ALGEBRAUTILS_H_
#include "Environment.h"
#include "RandomSeed.h"

#ifdef QPX
#include "AlgebraUtils_QPX.h"
#endif
#ifndef QPX

namespace Update {

class AlgebraUtils {
public:
	AlgebraUtils();
	~AlgebraUtils();

	/**
	 * This function returns back the squared norm of a dirac_vector
	 * \return the squared norm of vector
	 */
	template<typename dirac_vector_t> static long_real_t squaredNorm(const dirac_vector_t& vector) {
		long_real_t result = 0.;
#pragma omp parallel for reduction(+:result)
		for (int site = 0; site < vector.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				result += real(vector_dot(vector[site][mu],vector[site][mu]));
			}
		}
		reduceAllSum(result);
		return result;
	}

	/**
	 * This function normalizes the vector, norm(vector) = 1
	 */
	template<typename dirac_vector_t> static void normalize(dirac_vector_t& vector) {
		long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
		for (int site = 0; site < vector.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				norm += real(vector_dot(vector[site][mu],vector[site][mu]));
			}
		}
		reduceAllSum(norm);
#pragma omp parallel for
		for (int site = 0; site < vector.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				vector[site][mu] = vector[site][mu]/sqrt(norm);
			}
		}
	}

	/**
	 * This function generates a real random gaussian vector, normalized to exp(-x^2)
	 */
	template<typename dirac_vector_t> static void generateRandomGaussianVector(dirac_vector_t& vector) {
#ifndef MULTITHREADING
		random_generator_t rng(RandomSeed::randomSeed());
		random_normal_generator_t randomNormal(RandomSeed::getNormalNumberGenerator(rng,1./sqrt(2.)));
#endif
#ifdef MULTITHREADING
		random_generator_t** randomGenerator = new random_generator_t*[omp_get_max_threads()];
		random_normal_generator_t** randomNormal = new random_normal_generator_t*[omp_get_max_threads()];
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
			randomNormal[i] = new random_normal_generator_t(RandomSeed::getNormalNumberGenerator(*randomGenerator[i],1./sqrt(2.)));
		}
#endif
#pragma omp parallel for
		for (int site = 0; site < vector.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
					vector[site][mu][c] = std::complex<real_t>(randomNormal(), 0.);
#endif
#ifdef MULTITHREADING
					vector[site][mu][c] = std::complex<real_t>((*randomNormal[omp_get_thread_num()])(), 0.);
#endif
				}
			}
		}
		vector.updateHalo();
#ifdef MULTITHREADING
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			delete randomGenerator[i];
			delete randomNormal[i];
		}
		delete[] randomGenerator;
		delete[] randomNormal;
#endif
	}

	/**
	 * This function generates a random vector
	 */
	template<typename dirac_vector_t> static void generateRandomVector(dirac_vector_t& vector) {
#ifndef MULTITHREADING
		random_generator_t rng(RandomSeed::randomSeed());
		random_uniform_generator_t randomUniform(RandomSeed::getRandomNumberGenerator(rng));
#endif
#ifdef MULTITHREADING
		random_generator_t** randomGenerator = new random_generator_t*[omp_get_max_threads()];
		random_uniform_generator_t** randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
			randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
		}
#endif
#pragma omp parallel for
		for (int site = 0; site < vector.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
					vector[site][mu][c] = std::complex<real_t>(randomUniform(), randomUniform());
#endif
#ifdef MULTITHREADING
					vector[site][mu][c] = std::complex<real_t>((*randomUniform[omp_get_thread_num()])(), (*randomUniform[omp_get_thread_num()])());
#endif
				}
			}
		}
		vector.updateHalo();
#ifdef MULTITHREADING
		for (int i = 0; i < omp_get_max_threads(); ++i) {
			delete randomGenerator[i];
			delete randomUniform[i];
		}
		delete[] randomGenerator;
		delete[] randomUniform;
#endif
	}

	/**
	 * This function returns back the dot product of two dirac_vector
	 * \return vector1.vector.2
	 */
	template<typename dirac_vector_t> static std::complex<long_real_t> dot(const dirac_vector_t& vector1, const dirac_vector_t& vector2) {
		long_real_t result_re = 0.;
		long_real_t result_im = 0.;
#pragma omp parallel for reduction(+:result_re, result_im)
		for (int site = 0; site < vector1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(vector1[site][mu],vector2[site][mu]);
				result_re += real(partial);
				result_im += imag(partial);
			}
		}
		reduceAllSum(result_re);
		reduceAllSum(result_im);
		return std::complex<long_real_t>(result_re, result_im);
	}

	/**
	 * This function returns back the dot product of two dirac_vector, with the first multiplied by gamma5
	 * \return vector1.gamma5.vector.2
	 */
	template<typename dirac_vector_t> static std::complex<long_real_t> gamma5dot(const dirac_vector_t& vector1, const dirac_vector_t& vector2) {
		long_real_t result_re = 0.;
		long_real_t result_im = 0.;
#pragma omp parallel for reduction(+:result_re, result_im)
		for (int site = 0; site < vector1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 2; ++mu) {
				complex partial = vector_dot(vector1[site][mu],vector2[site][mu]);
				result_re += real(partial);
				result_im += imag(partial);
			}
			for (unsigned int mu = 2; mu < 4; ++mu) {
				complex partial = vector_dot(vector1[site][mu],vector2[site][mu]);
				result_re -= real(partial);
				result_im -= imag(partial);
			}
		}
		reduceAllSum(result_re);
		reduceAllSum(result_im);
		return std::complex<long_real_t>(result_re, result_im);
	}

	/**
	 * This function returns back the norm product of (vector1 - vector2)
	 * \return norm(vector1 - vector2)
	 */
	template<typename dirac_vector_t> static long_real_t differenceNorm(const dirac_vector_t& vector1, const dirac_vector_t& vector2) {
		long_real_t result = 0.;
#pragma omp parallel for reduction(+:result)
		for (int site = 0; site < vector1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				result += real(vector_dot((vector1[site][mu]-vector2[site][mu]),(vector1[site][mu]-vector2[site][mu])));
			}
		}
		reduceAllSum(result);
		return result;
	}

	template<typename dirac_vector_t> static void setToZero(dirac_vector_t& vector) {
#pragma omp parallel for
		for (int site = 0; site < vector.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_zero(vector[site][mu]);
			}
		}
	}
};

} /* namespace Update */

#endif

#endif /* ALGEBRAUTILS_H_ */
