#ifndef DIRACEIGENSOLVER_H_
#define DIRACEIGENSOLVER_H_
#include "dirac_operators/DiracOperator.h"
#include "inverters/BiConjugateGradient.h"
#include "dirac_functions/Polynomial.h"
#include "dirac_functions/ChebyshevRecursion.h"

namespace Update {

enum EigevaluesMode {LargestReal = 0, SmallestReal, LargestImaginary, SmallestImaginary};

class DiracEigenSolver {
	double epsilon;
	double inverterPrecision;
	unsigned int inverterMaximumSteps;
	unsigned int extra_steps;
	bool useChebyshev;
	unsigned int maximalNumberOfRestarts;
public:
	DiracEigenSolver();
	DiracEigenSolver(const DiracEigenSolver& copy);
	~DiracEigenSolver();

	void setInverterPrecision(const real_t& precision);
	real_t getInverterPrecision() const;

	void setTolerance(const real_t& precision);
	real_t getTolerance() const;

	void setExtraSteps(unsigned int _extra_steps);
	unsigned int getExtraSteps() const;

	void setMaximalNumberOfRestarts(unsigned int _maximalNumberOfRestarts);
	unsigned int getMaximalNumberOfRestarts() const;

	void setUseChebyshev(bool _useChebyshev);
	bool getUseChebyshev() const;

	void maximumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int n, EigevaluesMode mode = LargestReal);
	void minimumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors/*, Polynomial& map*/, unsigned int n/*, int nmode*/);
	void minimumNonHermitianEigenvalues(DiracOperator* diracOperator, DiracOperator* squareHermitianDiracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int n);

	void forceHermitianPairing(std::vector<reduced_dirac_vector_t>& V, std::vector< std::complex<real_t> >& eigenvalues);

	BiConjugateGradient* biConjugateGradient;

	ChebyshevRecursion* getChebyshevRecursion() const;
protected:
	ChebyshevRecursion* chebyshevRecursion;
	
	void extendArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, unsigned int m, unsigned int k, EigevaluesMode mode);
	void startArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, EigevaluesMode mode);
	long_real_t finishArnoldi(DiracOperator* diracOperator, const std::vector<reduced_dirac_vector_t>& V, const matrix_t& H, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int steps, EigevaluesMode mode);
	void restartArnoldi(DiracOperator* diracOperator, std::vector<reduced_dirac_vector_t>& V, reduced_dirac_vector_t& f,  matrix_t& H, unsigned int extra_steps);	
};

} /* namespace Update */
#endif /* DIRACEIGENSOLVER_H_ */
