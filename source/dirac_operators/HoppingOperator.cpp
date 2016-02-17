#include "HoppingOperator.h"

namespace  Update {

HoppingOperator::HoppingOperator(DiracOperator* _dirac) : dirac(_dirac) { }

void HoppingOperator::apply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input) {
	dirac->multiply(output, input);
	
#pragma omp parallel for
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			output[site][mu] = output[site][mu] - input[site][mu];
		}
	}
}

}
