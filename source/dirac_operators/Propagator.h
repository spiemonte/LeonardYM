/*
 * Propagator.h
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "LatticeSweep.h"
#include "DiracOperator.h"

namespace Update {

class Propagator {
public:
	static void constructPropagator(DiracOperator* diracOperator, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution);
};

}

#endif
