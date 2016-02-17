#include "MultiGridOperator.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

MultiGridOperator::MultiGridOperator() : dirac(0) { }

void MultiGridOperator::multiply(multigrid_vector_t& output, const multigrid_vector_t& input) {
	AlgebraUtils::setToZero(tmp_input);

#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif
	
	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) {
#pragma omp parallel for
		for (int site = 0; site < vectorspace[i]->localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					tmp_input[site][mu][c] += input(i,site)*(*vectorspace[i])[site][mu][c];
				}
			}
		}
	}

	tmp_input.updateHalo();

	dirac->multiply(tmp_output,tmp_input);//TODO: we don't need an updateHalo here

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
				proj[multigrid_vector_t::Layout::index(site)][processor] += vector_dot((*vectorspace[i])[site][mu],tmp_output[site][mu]);
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

void MultiGridOperator::multiplyAdd(multigrid_vector_t& output, const multigrid_vector_t& input, const complex& alpha) {
	AlgebraUtils::setToZero(tmp_input);

#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif
	
	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) {
#pragma omp parallel for
		for (int site = 0; site < vectorspace[i]->localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmp_input[site][mu] += input(i,site)*(*vectorspace[i])[site][mu];
			}
		}
	}
	
	tmp_input.updateHalo();
	dirac->multiplyAdd(tmp_output,tmp_input,tmp_input,alpha);

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
				proj[multigrid_vector_t::Layout::index(site)][processor] += vector_dot((*vectorspace[i])[site][mu],tmp_output[site][mu]);
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

matrix_t MultiGridOperator::asMatrix() {
	matrix_t result(multigrid_vector_t::Layout::size,multigrid_vector_t::Layout::size);
	multigrid_vector_t input, output;
	for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) {
#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			input[m] = 0.;
			if (m == i) input[m] = 1;
		}

		this->multiply(output,input);

#pragma omp parallel for
		for (int j = 0; j < multigrid_vector_t::Layout::size; ++j) {
			result(j,i) = output[j];
		}

	}

	return result;
}

}
