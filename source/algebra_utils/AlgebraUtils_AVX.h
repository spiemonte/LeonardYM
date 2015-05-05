#ifndef ALGEBRAUTILS_AVX_H_
#define ALGEBRAUTILS_AVX_H_
#include "Environment.h"
#include "utils/RandomSeed.h"

#ifdef AVX
//avxintrin.h
#include <immintrin.h>

typedef union {
			double value;
			char alignment [32];
		} aligned_double;

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
		/*int ds = vector1.localsize % 4;
		
#pragma omp parallel for reduction(+:result_re, result_im)
		for (int site = 0; site < vector1.localsize - ds; site += 4) {
			double tmp[] = {0.,0.,0.,0.};
			__m256d re = _mm256_load_pd(tmp);
			__m256d im = _mm256_load_pd(tmp);
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					double vector1real[4], vector2real[4], vector1imag[4], vector2imag[4];
					for (unsigned int j = 0; j < 4; ++j) {
						vector1real[j] = vector1[site+j][mu][c].real();
						vector1imag[j] = vector1[site+j][mu][c].imag();
						vector2real[j] = vector2[site+j][mu][c].real();
						vector2imag[j] = vector2[site+j][mu][c].imag();
					}
					__m256d v1r = _mm256_load_pd(vector1real);
					__m256d v1i = _mm256_load_pd(vector1imag);
					__m256d v2r = _mm256_load_pd(vector2real);
					__m256d v2i = _mm256_load_pd(vector2imag);
					re = _mm256_add_pd(re, _mm256_add_pd(_mm256_mul_pd(v1r,v2r),_mm256_mul_pd(v1i,v2i)));
					im = _mm256_add_pd(im, _mm256_sub_pd(_mm256_mul_pd(v1r,v2i),_mm256_mul_pd(v1i,v2r)));
				}
			}
			_mm256_store_pd(tmp,re);
			result_re += tmp[0] + tmp[1] + tmp[2] + tmp[3];
			_mm256_store_pd(tmp,im);
			result_im += tmp[0] + tmp[1] + tmp[2] + tmp[3];
		}


		for (int site = vector1.localsize - ds; site < vector1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(vector1[site][mu],vector2[site][mu]);
				result_re += real(partial);
				result_im += imag(partial);
			}
		}*/

		const double* ptr1 = reinterpret_cast<const double*>(vector1.getRawData());
		const double* ptr2 = reinterpret_cast<const double*>(vector2.getRawData());

		/*double vector1realr[] = {ptr1[0+0],ptr1[0+2],ptr1[0+4],ptr1[0+6]};
		double vector1imagr[] = {ptr1[0+1],ptr1[0+3],ptr1[0+5],ptr1[0+7]};

		std::cout << "vediamo: " << vector1realr[1] << " " << vector1imagr[1] << " " << vector1[0][0][1].real() << " " <<  vector1[0][0][1].imag() << std::endl;
		std::cout << "test: " << diracVectorLength*4*vector1.localsize << std::endl;*/

		//typedef double aligned_double __attribute__((aligned (32)));

#pragma omp parallel for reduction(+:result_re, result_im)
		for (int index = 0; index < 2*diracVectorLength*4*vector1.localsize; index += 16) {
			float vector1real[] __attribute__ ((aligned (16))) = {ptr1[index+0],ptr1[index+2],ptr1[index+4],ptr1[index+6],ptr1[index+8],ptr1[index+10],ptr1[index+12],ptr1[index+14]};
			float vector1imag[] __attribute__ ((aligned (16))) = {ptr1[index+1],ptr1[index+3],ptr1[index+5],ptr1[index+7],ptr1[index+9],ptr1[index+11],ptr1[index+13],ptr1[index+15]};
			float vector2real[] __attribute__ ((aligned (16))) = {ptr2[index+0],ptr2[index+2],ptr2[index+4],ptr2[index+6],ptr2[index+8],ptr2[index+10],ptr2[index+12],ptr2[index+14]};
			float vector2imag[] __attribute__ ((aligned (16))) = {ptr2[index+1],ptr2[index+3],ptr2[index+5],ptr2[index+7],ptr2[index+9],ptr2[index+11],ptr2[index+13],ptr2[index+15]};
			
			__m256 v1r = _mm256_load_ps(vector1real);
			__m256 v1i = _mm256_load_ps(vector1imag);
			__m256 v2r = _mm256_load_ps(vector2real);
			__m256 v2i = _mm256_load_ps(vector2imag);

			float tmp[8] __attribute__ ((aligned (16)));

			_mm256_store_ps(tmp, _mm256_add_ps(_mm256_mul_ps(v1r,v2r),_mm256_mul_ps(v1i,v2i)));
			result_re += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
			_mm256_store_ps(tmp, _mm256_sub_ps(_mm256_mul_ps(v1r,v2i),_mm256_mul_ps(v1i,v2r)));
			result_im += tmp[0] + tmp[1] + tmp[2] + tmp[3] + tmp[4] + tmp[5] + tmp[6] + tmp[7];
			/*result_re += ptr1[index]*ptr2[index] + ptr1[index+1]*ptr2[index+1];
			result_im += ptr1[index]*ptr2[index+1] - ptr1[index+1]*ptr2[index];*/
		}

		
		reduceAllSum(result_re);
		reduceAllSum(result_im);
		return std::complex<long_real_t>(result_re, result_im);
	}

	template<typename dirac_vector_t> static std::complex<long_real_t> slow_dot(const dirac_vector_t& vector1, const dirac_vector_t& vector2) {
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
		for (unsigned int site = 0; site < vector1.localsize; ++site) {
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
				vector[site][mu].zeros();
			}
		}
	}
};

} /* namespace Update */

#endif

#endif /* ALGEBRAUTILS_H_ */
