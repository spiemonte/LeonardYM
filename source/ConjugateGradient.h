/*
 * ConjugateGradient.h
 *
 *  Created on: Sep 27, 2012
 *      Author: spiem_01
 */

#ifndef CONJUGATEGRADIENT_H_
#define CONJUGATEGRADIENT_H_
#include "Environment.h"
#include "dirac_operators/DiracOperator.h"

namespace Update {

class ConjugateGradient {
public:
	ConjugateGradient();
	~ConjugateGradient();

	bool solvev(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

private:
	reduced_dirac_vector_t p;
	reduced_dirac_vector_t r;
	reduced_dirac_vector_t tmp;

	double epsilon;
	double lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;
};

} /* namespace Update */
#endif /* CONJUGATEGRADIENT_H_ */
