#ifndef ALGEBRAUTILS_H_
#define ALGEBRAUTILS_H_
#include "../Environment.h"
#include "utils/RandomSeed.h"

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
	 * This function generates a complex random gaussian vector, normalized to exp(-x^2)
	 */
	template<typename dirac_vector_t> static void generateRandomComplexGaussianVector(dirac_vector_t& vector) {
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
					vector[site][mu][c] = std::complex<real_t>(randomNormal(), randomNormal());
#endif
#ifdef MULTITHREADING
					vector[site][mu][c] = std::complex<real_t>((*randomNormal[omp_get_thread_num()])(), (*randomNormal[omp_get_thread_num()])());
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

	template<typename dirac_vector_t> static void conjugate(dirac_vector_t& output, const dirac_vector_t& input) {
#pragma omp parallel for
		for (int site = 0; site < input.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) output[site][mu][c] = conj(input[site][mu][c]);
			}
		}
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

	template<typename dirac_vector_t> static std::complex<long_real_t> real_dot(const dirac_vector_t& vector1, const dirac_vector_t& vector2) {
		long_real_t result_re = 0.;
		long_real_t result_im = 0.;
#pragma omp parallel for reduction(+:result_re, result_im)
		for (int site = 0; site < vector1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial(0.);
				for (unsigned int c = 0; c < diracVectorLength; ++c) partial += vector1[site][mu][c]*vector2[site][mu][c];
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

	/**
	 * This function multiplies the vector by gamma5
	 */
	template<typename dirac_vector_t> static void gamma5(dirac_vector_t& vector) {
#pragma omp parallel for
		for (int site = 0; site < vector.completesize; ++site) {
			for (unsigned int mu = 2; mu < 4; ++mu) {
				vector[site][mu] = - vector[site][mu];
			}
		}
	}

};

#ifdef ALIGNED_OPT

class AlignedAlgebraUtils {
public:
	template<typename aligned_vector_t> static void vector_plus_scalar_times_vector(
		aligned_vector_t& output,
		const aligned_vector_t& v1, 
		const std::complex<real_t>& alpha,
		const aligned_vector_t& v2) {

		real_t factor_real = alpha.real();
		real_t factor_imag = alpha.imag();

		const int localsize = v1.real_part[0].localsize;

#pragma vector aligned
#pragma omp parallel for simd
		for (int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) {
					const int i = 4*c + mu;
					output.real_part[i][site] = v1.real_part[i][site] + factor_real*v2.real_part[i][site] - factor_imag*v2.imag_part[i][site];
					output.imag_part[i][site] = v1.imag_part[i][site] + factor_real*v2.imag_part[i][site] + factor_imag*v2.real_part[i][site];
				}
			}
		}

	}

	template<typename aligned_vector_t> static real_t squaredNorm(const aligned_vector_t& v) {
		real_t result = 0.;
		const int localsize = v.real_part[0].localsize;

#pragma vector aligned
#pragma omp parallel for simd reduction(+:result)
		for (int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) {
					const int i = 4*c + mu;
					result += v.real_part[i][site]*v.real_part[i][site]+v.imag_part[i][site]*v.imag_part[i][site];
				}
			}
		}
		reduceAllSum(result);

		return result;
	}

	template<typename aligned_vector_t> static std::complex<real_t> dot(const aligned_vector_t& v1, const aligned_vector_t& v2) {
		real_t result_re = 0.;
		real_t result_im = 0.;
		const int localsize = v1.real_part[0].localsize;

#pragma vector aligned
#pragma omp parallel for simd reduction(+:result_re,result_im)
		for (int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) {
					const int i = 4*c + mu;
					result_re += v1.real_part[i][site]*v2.real_part[i][site]+v1.imag_part[i][site]*v2.imag_part[i][site];
					result_im += v1.real_part[i][site]*v2.imag_part[i][site]-v1.imag_part[i][site]*v2.real_part[i][site];
				}
			}
		}
		reduceAllSum(result_re);
		reduceAllSum(result_im);

		return std::complex<real_t>(result_re,result_im);
	}

	template<typename aligned_vector_t> static void setToZero(aligned_vector_t& v) {
		const int localsize = v.real_part[0].localsize;

#pragma vector aligned
#pragma omp parallel for simd
		for (int site = 0; site < localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) {
					const int i = 4*c + mu;
					v.real_part[i][site] = 0;
					v.imag_part[i][site] = 0;
				}
			}
		}
	}
};

#endif

} /* namespace Update */

#endif /* ALGEBRAUTILS_H_ */
