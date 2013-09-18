/*
 * ToString.h
 *
 *  Created on: May 24, 2012
 *      Author: spiem_01
 */

#ifndef TOSTRING_H_
#define TOSTRING_H_
#include "MatrixTypedef.h"
#include <sstream>
#include <iomanip>

namespace Update {

#ifdef EIGEN

template<class T, int r, int c> std::string toString(const Eigen::Matrix< T, r, c >& matrix) {
	std::ostringstream os;
	os << "{";
	if (c != 1) {
		for (unsigned int i = 0; i < r - 1; ++i) {
			os << "{";
			for (unsigned int j = 0; j < c - 1; ++j) {
				os << std::setprecision(11) << real(matrix(i,j)) << "+I*"<< imag(matrix(i,j)) << ",";
			}
			os << std::setprecision(11) << real(matrix(i,c-1)) << "+I*" << imag(matrix(i,c-1))  << "},";
		}
		os << std::setprecision(11) << "{";
		for (unsigned int j = 0; j < c - 1; ++j) {
			os << std::setprecision(11) << real(matrix(r-1,j)) << "+I*" << imag(matrix(r-1,j)) << ",";
		}
		os << std::setprecision(11) << real(matrix(r-1,c-1)) << "+I*" << imag(matrix(r-1,c-1)) << "}}";
	}
	else {
		for (unsigned int i = 0; i < r-1; ++i) {
			os << std::setprecision(11) << real(matrix(i,c-1)) << "+I*" << imag(matrix(i,c-1))  << ",";
		}
		os << std::setprecision(11) << real(matrix(r-1,c-1)) << "+I*" << imag(matrix(r-1,c-1)) << "}";
	}
	return os.str();
}

template<class T> std::string toString(const Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >& matrix) {
	std::ostringstream os;
	os << "{";
	unsigned int c = matrix.cols();
	unsigned int r = matrix.rows();
	for (unsigned int i = 0; i < r - 1; ++i) {
		os << "{";
		for (unsigned int j = 0; j < c - 1; ++j) {
			os << std::setprecision(11) << real(matrix(i,j)) << "+I*"<< imag(matrix(i,j)) << ",";
		}
		os << std::setprecision(11) << real(matrix(i,c-1)) << "+I*" << imag(matrix(i,c-1))  << "},";
	}
	os << std::setprecision(11) << "{";
	for (unsigned int j = 0; j < c - 1; ++j) {
		os << std::setprecision(11) << real(matrix(r-1,j)) << "+I*" << imag(matrix(r-1,j)) << ",";
	}
	os << std::setprecision(11) << real(matrix(r-1,c-1)) << "+I*" << imag(matrix(r-1,c-1)) << "}}";
	return os.str();
}

#endif

template<class T> std::string toString(const T& toConvert) {
	std::ostringstream os;
	os << toConvert;
	return os.str();
}

}



#endif /* TOSTRING_H_ */
