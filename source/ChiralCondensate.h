/*
 * ChiralCondensate.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef CHIRALCONDENSATE_H_
#define CHIRALCONDENSATE_H_

#include "LatticeSweep.h"
#include "StochasticEstimator.h"
#include "dirac_operators/DiracOperator.h"
#include "BiConjugateGradient.h"

namespace Update {

class ChiralCondensate : public Update::LatticeSweep, public Update::StochasticEstimator {
public:
	ChiralCondensate();
	ChiralCondensate(const ChiralCondensate& toCopy);
	~ChiralCondensate();

	virtual void execute(environment_t& environment);

private:
	extended_dirac_vector_t randomNoise;
	extended_dirac_vector_t tmp_square;
	extended_dirac_vector_t tmp;
	DiracOperator* squareDiracOperator;
	DiracOperator* diracOperator;
	BiConjugateGradient* biConjugateGradient;
};

} /* namespace Update */
#endif /* CHIRALCONDENSATE_H_ */
