#include "BlockBasis.h"
#include "MultiGridVector.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

BlockBasis::BlockBasis(int _basisDimension) : basisDimension(_basisDimension), vectorspace(new reduced_dirac_vector_t[_basisDimension]) {
	if (_basisDimension != multigrid_vector_t::Layout::basisDimension) {
		if (isOutputProcess()) std::cout << "BlockBasis::Reinitilization of the Layout " << std::endl;
		multigrid_vector_t::Layout::basisDimension = basisDimension;
		multigrid_vector_t::Layout::initialize();
	}
}

BlockBasis::BlockBasis(const BlockBasis& toCopy) : basisDimension(toCopy.basisDimension), vectorspace(new reduced_dirac_vector_t[toCopy.basisDimension]) {
	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) vectorspace[i] = toCopy.vectorspace[i];
}

BlockBasis::~BlockBasis() {
	delete[] vectorspace;
}

void BlockBasis::setBasisDimension(int _basisDimension) {
	if (_basisDimension != multigrid_vector_t::Layout::basisDimension) {
		multigrid_vector_t::Layout::basisDimension = basisDimension;
		multigrid_vector_t::Layout::initialize();
	}
	basisDimension = _basisDimension;
	delete[] vectorspace;
	vectorspace = new reduced_dirac_vector_t[_basisDimension];
}

void BlockBasis::orthogonalize(int index) {
#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif

	//We orthogonalize first the vector between the others to have a better projection to the low modes
	for (int j = 0; j < multigrid_vector_t::Layout::basisDimension; ++j) {
		if (j != index) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(vectorspace[j],vectorspace[index]));
#pragma omp parallel for
			for (int site = 0; site < vectorspace[index].completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					vectorspace[index][site][mu] -= proj*vectorspace[j][site][mu];
				}
			}
		}
	}
	//localBasis[i].updateHalo();TODO maybe not needed
	AlgebraUtils::normalize(vectorspace[index]);

	//Now we perform a block orthogonalization of the local basis
	for (int j = 0; j < multigrid_vector_t::Layout::basisDimension; ++j) {
		if (j != index) {
			std::complex<real_t> proj[multigrid_vector_t::Layout::totalNumberOfBlocks][numberProcessors];
			for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[k][p] = 0.;
			}
#pragma omp parallel for
			for (int site = 0; site < vectorspace[index].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[multigrid_vector_t::Layout::index(site)][0] += conj(vectorspace[j][site][mu][c])*vectorspace[index][site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[multigrid_vector_t::Layout::index(site)][omp_get_thread_num()] += conj(vectorspace[j][site][mu][c])*vectorspace[index][site][mu][c];
#endif
					}
				}
			}
			for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[k][0] += proj[k][p];
			}
			//for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[k][0]);
#pragma omp parallel for
			for (int site = 0; site < vectorspace[index].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					vectorspace[index][site][mu] -= proj[multigrid_vector_t::Layout::index(site)][0]*vectorspace[j][site][mu];
				}
			}
		}

		//Now we block normalize the vector
		real_t norm[multigrid_vector_t::Layout::totalNumberOfBlocks][numberProcessors];
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
			for (int p = 0; p < numberProcessors; ++p) norm[k][p] = 0.;
		}
#pragma omp parallel for
		for (int site = 0; site < vectorspace[index].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
					norm[multigrid_vector_t::Layout::index(site)][0] += real(conj(vectorspace[index][site][mu][c])*vectorspace[index][site][mu][c]);
#endif
#ifdef MULTITHREADING
					norm[multigrid_vector_t::Layout::index(site)][omp_get_thread_num()] += real(conj(vectorspace[index][site][mu][c])*vectorspace[index][site][mu][c]);
#endif
				}
			}
		}
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
			for (int p = 1; p < numberProcessors; ++p) norm[k][0] += norm[k][p];
		}
		//for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(norm[k][0]);
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) norm[k][0] = sqrt(norm[k][0]);
#pragma omp parallel for
		for (int site = 0; site < vectorspace[index].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				vectorspace[index][site][mu] = vectorspace[index][site][mu]/norm[multigrid_vector_t::Layout::index(site)][0];
			}
		}		

		//TODO: is it needed?
		vectorspace[index].updateHalo();
	}			
}

void BlockBasis::orthogonalize() {
#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif

	//We orthogonalize first the vectors between themself to have a better projection to the low modes
	for (int i = 0; i < basisDimension; ++i) {
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(vectorspace[j],vectorspace[i]));
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					vectorspace[i][site][mu] -= proj*vectorspace[j][site][mu];
				}
			}
		}
		//localBasis[i].updateHalo();TODO maybe not needed
		AlgebraUtils::normalize(vectorspace[i]);
	}

	//Now we perform a block orthogonalization of the local basis
	for (int i = 0; i < multigrid_vector_t::Layout::basisDimension; ++i) {
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj[multigrid_vector_t::Layout::totalNumberOfBlocks][numberProcessors];
			for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[k][p] = 0.;
			}
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[multigrid_vector_t::Layout::index(site)][0] += conj(vectorspace[j][site][mu][c])*vectorspace[i][site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[multigrid_vector_t::Layout::index(site)][omp_get_thread_num()] += conj(vectorspace[j][site][mu][c])*vectorspace[i][site][mu][c];
#endif
					}
				}
			}
			for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[k][0] += proj[k][p];
			}
			//for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[k][0]);
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					vectorspace[i][site][mu] -= proj[multigrid_vector_t::Layout::index(site)][0]*vectorspace[j][site][mu];
				}
			}
		}

		//Now we block normalize the vector
		real_t norm[multigrid_vector_t::Layout::totalNumberOfBlocks][numberProcessors];
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
			for (int p = 0; p < numberProcessors; ++p) norm[k][p] = 0.;
		}
#pragma omp parallel for
		for (int site = 0; site < vectorspace[i].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
					norm[multigrid_vector_t::Layout::index(site)][0] += real(conj(vectorspace[i][site][mu][c])*vectorspace[i][site][mu][c]);
#endif
#ifdef MULTITHREADING
					norm[multigrid_vector_t::Layout::index(site)][omp_get_thread_num()] += real(conj(vectorspace[i][site][mu][c])*vectorspace[i][site][mu][c]);
#endif
				}
			}
		}
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) {
			for (int p = 1; p < numberProcessors; ++p) norm[k][0] += norm[k][p];
		}
		//for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(norm[k][0]);
		for (int k = 0; k < multigrid_vector_t::Layout::totalNumberOfBlocks; ++k) norm[k][0] = sqrt(norm[k][0]);
#pragma omp parallel for
		for (int site = 0; site < vectorspace[i].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				vectorspace[i][site][mu] = vectorspace[i][site][mu]/norm[multigrid_vector_t::Layout::index(site)][0];
			}
		}		

		//TODO: is it needed?
		vectorspace[i].updateHalo();
	}
}

}

