#include "MultiGridProjector.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

MultiGridProjector::MultiGridProjector() { }

void MultiGridProjector::apply(multigrid_vector_t& output, const reduced_dirac_vector_t& input) {
#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif	
	
#pragma omp parallel for
	for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) {
		output[i] = 0.;
	}
	
	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) {
		//Initialize and set to zero the projections
		std::complex<real_t> proj[multigrid_vector_t::Layout::totalNumberOfBlocks][numberProcessors];
#pragma omp parallel for
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
			for (int p = 0; p < numberProcessors; ++p) proj[k][p] = 0.;
		}


#pragma omp parallel for
		for (int site = 0; site < vectorspace[i]->localsize; ++site) {
#ifndef MULTITHREADING
			int processor = 0;
#endif
#ifdef MULTITHREADING
			int processor = omp_get_thread_num();
#endif				
			for (unsigned int mu = 0; mu < 4; ++mu) {
				proj[multigrid_vector_t::Layout::index(site)][processor] += vector_dot((*vectorspace[i])[site][mu],input[site][mu]);
			}
		}

		for (int k = 0; k < numberProcessors; ++k) {
#pragma omp parallel for
			for (int j = 0; j < multigrid_vector_t::Layout::totalNumberOfBlocks; ++j) {
				output[i*multigrid_vector_t::Layout::totalNumberOfBlocks + j] += proj[j][k];
			}
		}
	}			
}

void MultiGridProjector::apply(reduced_dirac_vector_t& output, const multigrid_vector_t& input) {
	AlgebraUtils::setToZero(output);

	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) {
#pragma omp parallel for
		for (int site = 0; site < vectorspace[i]->localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += input(i,site)*(*vectorspace[i])[site][mu];
			}
		}
	}
	output.updateHalo();
}

}

