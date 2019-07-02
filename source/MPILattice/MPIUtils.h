#ifndef MPIUTILS_H
#define MPIUTILS_H

//Utils functions for MPI output to shell and for reduceAllSum

inline bool isOutputProcess() {
#ifdef ENABLE_MPI
	int this_processor;
	MPI_Comm_rank(MPI_COMM_WORLD, &this_processor);
	return this_processor == 0;
#endif
#ifndef ENABLE_MPI
	return true;
#endif
}

#ifdef ENABLE_MPI
inline void reduceAllSum(double& value) {
	double result = value;
	MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(double&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(float& value) {
	float result = value;
	MPI_Allreduce(&value, &result, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(float&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(int& value) {
	int result = value;
	MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(int&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(std::complex<double>& value) {
	double reValue = real(value), reResult;
	double imValue = imag(value), imResult;
	MPI_Allreduce(&reValue, &reResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&imValue, &imResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = std::complex<double>(reResult,imResult);
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(std::complex<double>&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(std::complex<long double>& value) {
        long double reValue = real(value), reResult;
        long double imValue = imag(value), imResult;
        MPI_Allreduce(&reValue, &reResult, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&imValue, &imResult, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        value = std::complex<long double>(reResult,imResult);
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(std::complex<long double>&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(std::complex<float>& value) {
	float reValue = real(value), reResult;
        float imValue = imag(value), imResult;
	MPI_Allreduce(&reValue, &reResult, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&imValue, &imResult, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	value = std::complex<float>(reResult,imResult);
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(std::complex<float>&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(long double& value) {
	long double result = value;
	MPI_Allreduce(&value, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(long double&) { }
#endif

#endif

