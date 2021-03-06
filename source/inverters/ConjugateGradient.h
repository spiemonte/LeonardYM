#ifndef CONJUGATEGRADIENT_H_
#define CONJUGATEGRADIENT_H_
#include "Environment.h"
#include "dirac_operators/DiracOperator.h"

namespace Update {

class ConjugateGradient {
public:
	ConjugateGradient();
	~ConjugateGradient();

	bool solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess = 0);

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
