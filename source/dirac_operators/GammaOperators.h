#ifndef GAMMA_OPERATORS_H
#define GAMMA_OPERATORS_H
#include "utils/Gamma.h"
#include "Environment.h"

namespace Update {

class GammaOperators {
public:
	
	void multiply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input, int index);
	
private:
	Gamma gammaMatrices;
};

}

#endif

