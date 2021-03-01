#ifndef MATRIXTYPEDEF_H_
#define MATRIXTYPEDEF_H_

//Matrix specification and typedefs for fundamental and adjoint reps

#define PI 3.1415926535897932384626433832795028841971693993751

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#ifdef MULTITHREADING
#include <omp.h>
#endif

#ifdef EIGEN

#define EIGEN_DISABLE_STACK_SIZE_ASSERT
#define NDEBUG
#define EIGEN_MATRIXBASE_PLUGIN MBE
#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <Eigen/QR>

namespace Update {

const int numberColors = NUMCOLORS;
typedef float single_real_t;
typedef double real_t;
typedef long double long_real_t;
typedef std::complex<real_t> complex;
typedef std::complex<single_real_t> single_complex;

inline std::complex<float> operator*(const double& x, const std::complex<float>& y) {
	return std::complex<float>(x*y.real(),x*y.imag());
}

inline std::complex<float> operator*(const std::complex<float>& y, const double& x) {
	return std::complex<float>(y.real()*x,y.imag()*x);
}

inline std::complex<float> operator/(const std::complex<float>& y, const double& x) {
	return std::complex<float>(y.real()/x,y.imag()/x);
}

inline std::complex<float> operator+(const std::complex<float>& y, const double& x) {
	return std::complex<float>(y.real()+x,y.imag());
}

inline std::complex<float> operator+(const double& x, const std::complex<float>& y) {
	return std::complex<float>(x+y.real(),y.imag());
}

typedef Eigen::Matrix< complex, 2, 2 > matrix2x2_t;
typedef Eigen::Matrix< complex, Eigen::Dynamic, Eigen::Dynamic > matrix_t;
typedef Eigen::Matrix< complex, Eigen::Dynamic, 1 > vector_t;

typedef Eigen::Matrix< single_complex, 2, 2 > single_matrix2x2_t;
typedef Eigen::Matrix< single_complex, Eigen::Dynamic, Eigen::Dynamic > single_matrix_t;

typedef Eigen::Matrix< complex, numberColors, numberColors > FundamentalGroup;
typedef Eigen::Matrix< real_t, numberColors*numberColors - 1, numberColors*numberColors - 1 > AdjointGroup;

typedef Eigen::Matrix< complex, numberColors, 1 > FundamentalVector;
typedef Eigen::Matrix< complex, numberColors*numberColors - 1, 1 > AdjointComplexVector;
typedef Eigen::Matrix< real_t, numberColors*numberColors - 1, 1 > AdjointRealVector;

typedef Eigen::Matrix< single_complex, numberColors, 1 > single_FundamentalVector;
typedef Eigen::Matrix< single_complex, numberColors*numberColors - 1, 1 > single_AdjointVector;
typedef Eigen::Matrix< real_t, numberColors*numberColors - 1, 1 > ForceVector;

#ifdef ADJOINT
typedef Eigen::Matrix< complex, numberColors, numberColors > GaugeGroup;
typedef Eigen::Matrix< real_t, numberColors*numberColors - 1, numberColors*numberColors - 1 > FermionicGroup;
typedef Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 > FermionicForceMatrix;
typedef Eigen::Matrix< complex, numberColors*numberColors - 1, 1 > GaugeVector;
typedef Eigen::Matrix< single_complex, numberColors*numberColors - 1, 1 > single_GaugeVector;

const int diracVectorLength = numberColors*numberColors - 1;
#endif
#ifndef ADJOINT
typedef Eigen::Matrix< complex, numberColors, numberColors > GaugeGroup;
typedef Eigen::Matrix< complex, numberColors, numberColors > FermionicGroup;
typedef Eigen::Matrix< complex, numberColors, numberColors > FermionicForceMatrix;
typedef Eigen::Matrix< complex, numberColors, 1 > GaugeVector;
typedef Eigen::Matrix< single_complex, numberColors, 1 > single_GaugeVector;

const int diracVectorLength = numberColors;
#endif

typedef Eigen::Matrix< complex, 2, 2 > FundamentalSU2Group;
typedef Eigen::Matrix< real_t, 3, 3 > AdjointSU2Group;

#define htrans(x) ((x).adjoint())
#define dag(x) ((x).conjugate())
#define adj(x) ((x).adjoint())
#define det(x) (x).determinant()
#define trace(x) (x).trace()
#define inverse(x) (x).inverse()
#define n_rows rows()
#define n_cols cols()
#define set_to_zero(x) ((x).zeros())
#define set_to_identity(x) ((x) = GaugeGroup::Identity())
#define vector_dot(x,y) (x).dot(y)

#define qr(q,r,x) Eigen::HouseholderQR<GaugeGroup> qrdecomp(x); \
q = qrdecomp.householderQ()

}

#endif
#ifdef ARMADILLO

#define ARMA_NO_DEBUG
//#define ARMA_DONT_USE_BLAS
#include <armadillo>

namespace Update {

const int numberColors = NUMCOLORS;
typedef double real_t;
typedef double long_real_t;
typedef std::complex<real_t> complex;

typedef arma::Mat< complex >::fixed<2,2> matrix2x2_t;
typedef arma::Mat< complex > matrix_t;
typedef arma::Col< complex > vector_t;

typedef arma::Mat< complex >::fixed<numberColors,numberColors> FundamentalGroup;
typedef arma::Mat< real_t >::fixed<numberColors*numberColors - 1, numberColors*numberColors - 1 > AdjointGroup;

typedef arma::Col< complex >::fixed< numberColors > FundamentalVector;
typedef arma::Col< complex >::fixed< numberColors*numberColors - 1 > AdjointVector;

#ifdef ADJOINT
typedef arma::Mat< complex >::fixed<numberColors,numberColors> GaugeGroup;
typedef arma::Mat< real_t >::fixed<numberColors*numberColors - 1, numberColors*numberColors - 1 > FermionicGroup;
typedef arma::Mat< complex >::fixed<numberColors,numberColors> FermionicForceMatrix;
typedef arma::Col< complex >::fixed< numberColors*numberColors - 1 > GaugeVector;

const int diracVectorLength = numberColors*numberColors - 1;
#endif
#ifndef ADJOINT
typedef arma::Mat< complex >::fixed<numberColors,numberColors> GaugeGroup;
typedef arma::Mat< complex >::fixed<numberColors,numberColors> FermionicGroup;
typedef arma::Mat< complex >::fixed<numberColors,numberColors> FermionicForceMatrix;
typedef arma::Col< complex >::fixed< numberColors > GaugeVector;

const int diracVectorLength = numberColors;
#endif

typedef arma::Mat< complex >::fixed<2,2> FundamentalSU2Group;
typedef arma::Mat< real_t >::fixed<3,3> AdjointSU2Group;

#define set_to_zero(x) ((x).zeros())
#define set_to_identity(x) (x).zeros(); (x.diag() += 1)
#define vector_dot(x,y) cdot((x),(y))
#define inverse(x) (x).i()
#define adj(x) (htrans(x))

}

#endif
#ifdef MATRIXTOOLKIT

#include "utils/MatrixToolkit.h"

namespace Update {

const int numberColors = NUMCOLORS;
typedef double real_t;
typedef long double long_real_t;
typedef std::complex<real_t> complex;
typedef std::complex<float> single_complex;

typedef matrix_toolkit::Matrix< complex, 2 > matrix2x2_t;
typedef matrix_toolkit::Matrix< complex, -1 > matrix_t;
typedef matrix_toolkit::Vector< complex, -1 > vector_t;

typedef matrix_toolkit::Matrix< complex, numberColors > FundamentalGroup;
typedef matrix_toolkit::Matrix< real_t, numberColors*numberColors - 1 > AdjointGroup;

typedef matrix_toolkit::Vector< complex, numberColors > FundamentalVector;
typedef matrix_toolkit::Vector< complex, numberColors*numberColors - 1 > AdjointVector;

#ifdef ADJOINT
typedef matrix_toolkit::Matrix< complex, numberColors > GaugeGroup;
typedef matrix_toolkit::Matrix< real_t, numberColors*numberColors - 1 > FermionicGroup;
typedef matrix_toolkit::Matrix< complex, numberColors*numberColors - 1 > FermionicForceMatrix;
typedef matrix_toolkit::Vector< complex, numberColors*numberColors - 1 > GaugeVector;
typedef matrix_toolkit::Vector< single_complex, numberColors*numberColors - 1 > single_GaugeVector;

const int diracVectorLength = numberColors*numberColors - 1;
#endif
#ifndef ADJOINT
typedef matrix_toolkit::Matrix< complex, numberColors > GaugeGroup;
typedef matrix_toolkit::Matrix< complex, numberColors > FermionicGroup;
typedef matrix_toolkit::Matrix< complex, numberColors > FermionicForceMatrix;
typedef matrix_toolkit::Vector< complex, numberColors > GaugeVector;
typedef matrix_toolkit::Vector< single_complex, numberColors, 1 > single_GaugeVector;

const int diracVectorLength = numberColors;
#endif

typedef matrix_toolkit::Matrix< complex, 2 > FundamentalSU2Group;
typedef matrix_toolkit::Matrix< real_t, 3 > AdjointSU2Group;

}

#endif

#endif /* MATRIXTYPEDEF_H_ */
