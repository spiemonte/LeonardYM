#include "Propagator.h"
#include "OverlapOperator.h"

namespace Update {

void Propagator::constructPropagator(DiracOperator* diracOperator, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution) {
	if (diracOperator->getName() == "Overlap") {
		OverlapOperator* overlap = dynamic_cast<OverlapOperator*>(diracOperator);	

		real_t mass = overlap->getMass();
		overlap->setMass(0);

		reduced_dirac_vector_t tmp;
		
		overlap->multiply(tmp,solution);

#pragma omp parallel for
		for (int site = 0; site < tmp.completesize; ++site) {
			for (unsigned int mu = 0; mu< 4; ++mu) {
				solution[site][mu] = solution[site][mu] - tmp[site][mu];
			}
		}

		/*reduced_dirac_vector_t mah = tmp;
		AlgebraUtils::gamma5(mah);

		AlgebraUtils::gamma5(randomNoise);

		inverter->solve(diracOperator, randomNoise, tmp);

		if (diracOperator->getName() == "Overlap") {
#pragma omp parallel for
			for (int site = 0; site < tmp.completesize; ++site) {
				for (unsigned int mu = 0; mu< 4; ++mu) {
					tmp[site][mu] = tmp[site][mu] - randomNoise[site][mu];
				}
			}
#pragma omp parallel for
			for (int site = 0; site < tmp.completesize; ++site) {
				for (unsigned int mu = 0; mu< 4; ++mu) {
					tmp[site][mu] = tmp[site][mu] + mah[site][mu];
				}
			}
			long_real_t test = AlgebraUtils::squaredNorm(tmp);

			std::cout << "Giusto per: " << test << std::endl;
		}*/

		overlap->setMass(mass);
	}
}

}

