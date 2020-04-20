#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include "LatticeSweep.h"
#include "DiracOperator.h"

namespace Update {

class Propagator {
public:
#ifdef ENABLE_MPI
	static void constructPropagator(DiracOperator* diracOperator, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution);
#endif

	static void constructPropagator(DiracOperator* diracOperator, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution);
};

}

#endif
