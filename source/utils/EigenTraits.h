/*
 * EigenTraits.h
 *
 *  Created on: Mar 16, 2012
 *      Author: spiem_01
 */

#ifndef EIGENTRAITS_H_
#define EIGENTRAITS_H_
//#include "../math/grouptraits.h"
#include "MatrixTypedef.h"

namespace Update {

template<class T> class Group_Traits {
};

template<> class Group_Traits<Update::GaugeGroup> {
public:
	typedef Update::real_t Real;
	typedef Update::GaugeGroup Matrix;
	typedef Update::complex BasicElement;
	//	typedef Update::GaugeGroup MatrixForm;
	typedef Update::GaugeGroup BasicForm;
	typedef Update::GaugeGroup& Adaptor;
	typedef const Update::GaugeGroup& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors, 1 > RealVector;
	typedef Real RealVectorElement;
	typedef Eigen::Matrix< Update::complex, Update::numberColors, 1 > ComplexVector;
	typedef BasicElement ComplexVectorElement;
	static const size_t vectorsize = Update::numberColors*Update::numberColors;
	static const size_t matrixsize = Update::numberColors;
	static const size_t matrixdatalength = matrixsize*matrixsize;
	static const size_t vectordatalength = vectorsize;
	static const size_t matrixbytesize = matrixdatalength*sizeof(BasicElement);
	static const size_t bytesize = matrixdatalength*sizeof(BasicElement);
	typedef Eigen::Matrix< BasicElement, Update::numberColors, 1 > BasicVector;
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}
};


template<> class Group_Traits<Update::AdjointGroup> {
public:
	typedef Update::real_t Real;
	typedef Update::AdjointGroup Matrix;
	typedef Update::real_t BasicElement;
	typedef Update::AdjointGroup MatrixForm;
	typedef Update::AdjointGroup BasicForm;
	typedef Update::AdjointGroup& Adaptor;
	typedef const Update::AdjointGroup& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors*Update::numberColors - 1, 1 > RealVector;
	typedef Eigen::Matrix< Update::complex, Update::numberColors*Update::numberColors - 1, 1 > ComplexVector;
	static const size_t vectorsize = (Update::numberColors*Update::numberColors - 1)*(Update::numberColors*Update::numberColors - 1);
	static const size_t matrixsize = (Update::numberColors*Update::numberColors - 1)*(Update::numberColors*Update::numberColors - 1);
	static const size_t datalength = (Update::numberColors*Update::numberColors - 1)*(Update::numberColors*Update::numberColors - 1);//TODO
	static const size_t bytesize = datalength * sizeof(BasicElement);
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}
};

#if NUMCOLORS > 2
template<> class Group_Traits<Update::FundamentalSU2Group> {
public:
	typedef Update::real_t Real;
	typedef Update::GaugeGroup Matrix;
	typedef Update::complex BasicElement;
	typedef Update::GaugeGroup MatrixForm;
	typedef Update::GaugeGroup BasicForm;
	typedef Update::GaugeGroup& Adaptor;
	typedef const Update::GaugeGroup& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors, 1 > RealVector;
	typedef Eigen::Matrix< Update::complex, Update::numberColors, 1 > ComplexVector;
	static const size_t matrixsize = Update::numberColors*Update::numberColors;
	static const size_t datalength = Update::numberColors*Update::numberColors;//TODO
	static const size_t bytesize = datalength * sizeof(BasicElement);
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}
};


template<> class Group_Traits<Update::AdjointSU2Group> {
public:
	typedef Update::real_t Real;
	typedef Update::GaugeGroup Matrix;
	typedef Update::real_t BasicElement;
	typedef Update::GaugeGroup MatrixForm;
	typedef Update::GaugeGroup BasicForm;
	typedef Update::GaugeGroup& Adaptor;
	typedef const Update::GaugeGroup& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors, 1 > RealVector;
	typedef Eigen::Matrix< Update::complex, Update::numberColors, 1 > ComplexVector;
	static const size_t matrixsize = Update::numberColors*Update::numberColors;
	static const size_t datalength = Update::numberColors*Update::numberColors;//TODO
	static const size_t bytesize = datalength * sizeof(BasicElement);
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}
};

#endif


//TODO errato per adjoint
template<> class Group_Traits<Update::FundamentalVector> {
public:
	typedef Update::real_t Real;
	typedef Update::GaugeVector Matrix;
	typedef Update::complex BasicElement;
	typedef Update::GaugeVector MatrixForm;
	typedef Update::GaugeVector BasicForm;
	typedef Update::GaugeVector& Adaptor;
	typedef const Update::GaugeVector& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors, 1 > RealVector;
	typedef Eigen::Matrix< Update::complex, Update::numberColors, 1 > ComplexVector;
	static const size_t vectorsize = Update::numberColors;
	static const size_t matrixsize = Update::numberColors;
	static const size_t datalength = Update::numberColors;//TODO
	static const size_t bytesize = datalength * sizeof(BasicElement);
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}

	typedef Update::complex ComplexVectorElement;
};

template<> class Group_Traits<Update::AdjointVector> {
public:
	typedef Update::real_t Real;
	typedef Update::GaugeVector Matrix;
	typedef Update::complex BasicElement;
	typedef Update::GaugeVector MatrixForm;
	typedef Update::GaugeVector BasicForm;
	typedef Update::GaugeVector& Adaptor;
	typedef const Update::GaugeVector& const_Adaptor;
	typedef Eigen::Matrix< Update::real_t, Update::numberColors*Update::numberColors - 1, 1 > RealVector;
	typedef Eigen::Matrix< Update::complex, Update::numberColors*Update::numberColors - 1, 1 > ComplexVector;
	static const size_t vectorsize = Update::numberColors*Update::numberColors - 1;
	static const size_t matrixsize = Update::numberColors*Update::numberColors - 1;
	static const size_t datalength = Update::numberColors*Update::numberColors - 1;//TODO
	static const size_t bytesize = sizeof(Update::AdjointVector);//datalength * sizeof(BasicElement);
	static void setToZero(Matrix& m) {
		m = m.Zero();
	}
	static void setToIdentity(Matrix& m) {
		m = m.Identity();
	}

	typedef Update::complex ComplexVectorElement;
};

}

#endif /* EIGENTRAITS_H_ */
