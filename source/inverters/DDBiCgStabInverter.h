/*
 * DDBiCgStabInverter.h
 *
 *  Created on: Feb 26, 2013
 *      Author: spiem_01
 */

#ifndef DDBICGSTABINVERTER_H_
#define DDBICGSTABINVERTER_H_
#include "../dirac_operators/DiracOperator.h"

namespace Update {

class DDBiCgStabInverter {
public:
	DDBiCgStabInverter(unsigned int numberLowModeProjection);
	~DDBiCgStabInverter();

	bool solve(DiracOperator* dirac, const dirac_vector_t& source, dirac_vector_t& solution);

	void setPrecision(double _epsilon);
	double getPrecision() const;

	double getLastError() const;
	unsigned int getLastSteps() const;

	void setMaximumSteps(unsigned int _maxSteps);
	unsigned int getMaximumSteps() const;

private:
	std::vector< dirac_vector_t > low_modes;

	double epsilon;
	double lastError;
	unsigned int lastSteps;
	unsigned int maxSteps;

	dirac_vector_t residual;
	dirac_vector_t residual_hat;
	dirac_vector_t p;
	dirac_vector_t nu;
	dirac_vector_t s;
	dirac_vector_t t;

};

} /* namespace Update */
#endif /* DDBICGSTABINVERTER_H_ */
