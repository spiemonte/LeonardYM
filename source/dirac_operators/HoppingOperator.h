#ifndef HOPPINGOPERATOR_H
#define HOPPINGOPERATOR_H
#include "DiracOperator.h"

namespace Update {

class HoppingOperator {
public:
	HoppingOperator(DiracOperator* _dirac);

	void apply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input);
	
private:
	DiracOperator* dirac;
};

}

#endif

