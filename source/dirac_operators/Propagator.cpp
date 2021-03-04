#include "Propagator.h"
#include "OverlapOperator.h"

namespace Update {

#ifdef ENABLE_MPI
void Propagator::constructPropagator(DiracOperator* diracOperator, const extended_dirac_vector_t& source, extended_dirac_vector_t& solution) {
        if (diracOperator->getName() == "Overlap" || diracOperator->getName() == "ExactOverlap") {
                OverlapOperator* overlap = dynamic_cast<OverlapOperator*>(diracOperator);

                real_t mass = overlap->getMass();
                overlap->setMass(0.);

                extended_dirac_vector_t tmp;

                diracOperator->multiply(tmp,source);

#pragma omp parallel for
                for (int site = 0; site < tmp.completesize; ++site) {
                        for (unsigned int mu = 0; mu< 4; ++mu) {
                                solution[site][mu] = source[site][mu] - tmp[site][mu];
                        }
                }

		overlap->setMass(mass);
	}
	else {
		solution = source;
	}
}
#endif

void Propagator::constructPropagator(DiracOperator* diracOperator, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution) {
	if (diracOperator->getName() == "Overlap" || diracOperator->getName() == "ExactOverlap") {
		//Only for the overlap operator, we need to add a (1 - D_ov) factor
		OverlapOperator* overlap = dynamic_cast<OverlapOperator*>(diracOperator);	

		real_t mass = overlap->getMass();
		overlap->setMass(0.);

		reduced_dirac_vector_t tmp;
		
		overlap->multiply(tmp,source);

#pragma omp parallel for
		for (int site = 0; site < tmp.completesize; ++site) {
			for (unsigned int mu = 0; mu< 4; ++mu) {
				solution[site][mu] = source[site][mu] - tmp[site][mu];
			}
		}

		overlap->setMass(mass);
	}
	else {
		solution = source;
	}
}

}

