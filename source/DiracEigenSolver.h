/*
 * DiracEigenSolver.h
 *
 *  Created on: Jun 26, 2012
 *      Author: spiem_01
 */

#ifndef DIRACEIGENSOLVER_H_
#define DIRACEIGENSOLVER_H_
#include "dirac_operators/DiracOperator.h"
#include "BiConjugateGradient.h"
#include "Polynomial.h"

namespace Update {

class DiracEigenSolver {
	double epsilon;
	unsigned int extra_steps;
public:
	DiracEigenSolver();
	~DiracEigenSolver();

	void setPrecision(const real_t& precision);
	real_t getPrecision() const;

	void setExtraSteps(unsigned int _extra_steps);
	unsigned int getExtraSteps() const;

	void maximumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, unsigned int n);
	void minimumEigenvalues(DiracOperator* diracOperator, std::vector< std::complex<real_t> >& eigenvalues, std::vector<reduced_dirac_vector_t>& eigenvectors, Polynomial& map, unsigned int n, int nmode);

	BiConjugateGradient* biConjugateGradient;
};

} /* namespace Update */
#endif /* DIRACEIGENSOLVER_H_ */
