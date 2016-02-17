#include "GammaOperators.h"

namespace Update {

void GammaOperators::multiply(extended_dirac_vector_t& output, const extended_dirac_vector_t& input, int index) {
#pragma omp parallel for
	for (int site = 0; site < input.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(output[site][mu]);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				output[site][mu] += gammaMatrices.gammaChromaMatrices(index).at(mu,nu)*input[site][nu];
			}
		}
	}
}

}
