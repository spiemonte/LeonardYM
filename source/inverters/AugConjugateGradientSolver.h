/*
 * AugConjugateGradientSolver.h
 *
 *  Created on: Nov 4, 2013
 *      Author: spiem_01
 */

#ifndef AUGCONJUGATEGRADIENT_H_
#define AUGCONJUGATEGRADIENT_H_

namespace Update {

class AugConjugateGradient {
public:
	AugConjugateGradient();
	~AugConjugateGradient();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

private:
	double epsilon;
	double lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;

	reduced_dirac_vector_t residual;
	reduced_dirac_vector_t zeta;
	reduced_dirac_vector_t p;
	reduced_dirac_vector_t tmp;
	int m;
	std::complex<real_t>* wproj;
	reduced_dirac_vector_t* w;
	reduced_dirac_vector_t* dw;
};

} /* namespace Update */
#endif /* AUGCONJUGATEGRADIENTSOLVER_H_ */
