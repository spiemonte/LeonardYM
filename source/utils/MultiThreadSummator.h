#ifndef MULTITHREADSUMMATOR_H
#define MULTITHREADSUMMATOR_H

namespace Update {

template<typename T> class MultiThreadSummator {
	public:
		MultiThreadSummator() : result(0.0) {
#ifdef ENABLE_OMP
			data = new T[omp_get_max_threads()];
			for (int i = 0; i < omp_get_max_threads(); ++i) data[i] = 0;
#else
			data = 0;
#endif
		}

		MultiThreadSummator(const MultiThreadSummator& toCopy) {
#ifdef ENABLE_OMP
			data = new T[omp_get_max_threads()];
			for (int i = 0; i < omp_get_max_threads(); ++i) data[i] = toCopy.data[i];
#else
			data = toCopy.data;
#endif
			result = toCopy.result;
		}

		~MultiThreadSummator() {
#ifdef ENABLE_OMP
			delete[] data;
#endif
		}

		inline void add(const T& value) {
#ifdef ENABLE_OMP
			data[omp_get_thread_num()] += value;
#else
			data += value;
#endif
		}

		

		void reset() {
#ifdef ENABLE_OMP
			for (int i = 0; i < omp_get_max_threads(); ++i) data[i] = 0;
#else
			data = 0;
#endif
			result = 0.0;
		}

		//Compute and return the result of the sum
		T computeResult() {
#ifdef ENABLE_OMP
			result = T(0.0);
			for (int i = 0; i < omp_get_max_threads(); ++i) result += data[i];

#ifdef ENABLE_MPI
			this->reduceAllSum(result);
			
#endif
			return result;
#else
#ifdef ENABLE_MPI
			result = data;
			this->reduceAllSum(result);
			return result;
#else
			return data;
#endif
#endif
		}

		//Return only the result stored inside the class previously computed by computeResult()
		T getResult() const {
			return result;
		}
		
	private:
#ifdef ENABLE_OMP
		T* data;
#else
		T data;
#endif
		T result;

	protected:
#ifdef ENABLE_MPI
		void reduceAllSum(std::complex<double>& value) const {
			double valueRe = value.real(), resultRe = value.real();
			double valueIm = value.imag(), resultIm = value.imag();
			MPI_Allreduce(&valueRe, &resultRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&valueIm, &resultIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			value = std::complex<double>(resultRe,resultIm);
		}

		void reduceAllSum(std::complex<float>& value) const {
			float valueRe = value.real(), resultRe = value.real();
                        float valueIm = value.imag(), resultIm = value.imag();
			MPI_Allreduce(&valueRe, &resultRe, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(&valueIm, &resultIm, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
			value = std::complex<float>(resultRe,resultIm);
		}

		void reduceAllSum(long double& value) const {
			long double result = value;
			MPI_Allreduce(&value, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			value = result;
		}

		void reduceAllSum(double& value) const {
			double result = value;
			MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			value = result;
		}
	
		void reduceAllSum(float& value) const {
			double result = value;
			MPI_Allreduce(&value, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
			value = result;
		}
		
#endif
};

}

#endif
