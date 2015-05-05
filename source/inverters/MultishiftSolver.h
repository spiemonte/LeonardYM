/*
 * MultishiftSolver.h
 *
 *  Created on: May 11, 2012
 *      Author: spiem_01
 */

#ifndef MULTISHIFTSOLVER_H_
#define MULTISHIFTSOLVER_H_
#include "dirac_operators/DiracOperator.h"

namespace Update {

class MultishiftSolver {
public:
	MultishiftSolver(real_t _epsilon, unsigned int _maxSteps);
	virtual ~MultishiftSolver();

	static MultishiftSolver* getInstance(const std::string& name);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	void setMaxSteps(unsigned int _maxSteps);
	unsigned int getMaxSteps() const;

	/**
	 * This function implements the multishift solver for the operator dirac (NB: it must hermitian and definite positive)
	 * @param dirac the dirac operator
	 * @param source
	 * @param solutions
	 * @param shifts
	 * @return false if the solver fails
	 */
	virtual bool solve(DiracOperator* dirac, const extended_dirac_vector_t& source, std::vector<extended_dirac_vector_t>& solutions, const std::vector<real_t>& shifts) = 0;

protected:
	double epsilon;
	unsigned int maxSteps;
};

} /* namespace Update */
#endif /* MULTISHIFTSOLVER_H_ */
