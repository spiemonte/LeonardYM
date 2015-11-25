/*
 * GMRESR.h
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#ifndef GMRESR_H_
#define GMRESR_H_
#include "dirac_operators/DiracOperator.h"
#include <vector>

namespace Update {

class GMRESR {
	double epsilon;
	double lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;
	unsigned int conjugateSpaceDimension;

	reduced_dirac_vector_t* c;
	reduced_dirac_vector_t* u;
	reduced_dirac_vector_t r;
	
public:
	GMRESR();
	GMRESR(const GMRESR& toCopy);
	~GMRESR();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner = 0, reduced_dirac_vector_t const* initial_guess = 0);

#ifdef ENABLE_MPI
	bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution, DiracOperator* preconditioner = 0, extended_dirac_vector_t const* initial_guess = 0);
#endif

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

	void setConjugateSpaceDimension(unsigned int _maxSteps);
	unsigned int getConjugateSpaceDimension() const;
};

} /* namespace Update */
#endif /* BICONJUGATEGRADIENT_H_ */
