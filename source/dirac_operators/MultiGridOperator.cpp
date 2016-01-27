#include "MultiGridOperator.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

int MultiGridVectorLayout::xBlockSize = 4;
int MultiGridVectorLayout::yBlockSize = 4;
int MultiGridVectorLayout::zBlockSize = 4;
int MultiGridVectorLayout::tBlockSize = 4;

int MultiGridVectorLayout::totalNumberOfBlocks;
int MultiGridVectorLayout::basisDimension = 40;
int MultiGridVectorLayout::size;
	
reduced_index_lattice_t* MultiGridVectorLayout::blockIndex;

void MultiGridVectorLayout::initialize() {
	typedef reduced_index_lattice_t::Layout LT;
	//number of blocks in the x,y,z,t direction
	int numberBX = LT::loc_x/xBlockSize + ((LT::loc_x % xBlockSize) != 0 ? 1 : 0);
	int numberBY = LT::loc_y/yBlockSize + ((LT::loc_y % yBlockSize) != 0 ? 1 : 0);
	int numberBZ = LT::loc_z/zBlockSize + ((LT::loc_z % zBlockSize) != 0 ? 1 : 0);
	int numberBT = LT::loc_t/tBlockSize + ((LT::loc_t % tBlockSize) != 0 ? 1 : 0);

	if (LT::loc_x % xBlockSize != 0 || LT::loc_y % yBlockSize != 0 || LT::loc_z % zBlockSize != 0 || LT::loc_t % tBlockSize != 0) {
		if (isOutputProcess()) std::cout << "MultiGridOperator::Warning, block grid does not evenly match the processor grid: (" << LT::loc_x % xBlockSize << "," << LT::loc_y % yBlockSize << "," << LT::loc_z % zBlockSize << "," << LT::loc_t % tBlockSize << ")" << std::endl;
	}

	totalNumberOfBlocks = numberBX*numberBY*numberBZ*numberBT;
	size = totalNumberOfBlocks*basisDimension;

	blockIndex = new reduced_index_lattice_t;
	
	for (int site = 0; site < blockIndex->localsize; ++site) {
		int x = (LT::globalIndexX(site) % LT::loc_x)/xBlockSize;
		int y = (LT::globalIndexY(site) % LT::loc_y)/yBlockSize;
		int z = (LT::globalIndexZ(site) % LT::loc_z)/zBlockSize;
		int t = (LT::globalIndexT(site) % LT::loc_t)/tBlockSize;

		(*blockIndex)[site] = numberBT*(numberBZ*(numberBY*x + y) + z) + t;
	}
	blockIndex->updateHalo();
}

void MultiGridVectorLayout::setBasisDimension(int _basisDimension) {
	basisDimension = _basisDimension;
	size = basisDimension*totalNumberOfBlocks;
}





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





BlockBasis::BlockBasis(unsigned int _basisDimension) : basisDimension(_basisDimension), vectorspace(new reduced_dirac_vector_t[_basisDimension]) {
	if (basisDimension != multigrid_vector_t::Layout::basisDimension) {
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

void BlockBasis::orthogonalize(unsigned int index) {
#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif

	//We orthogonalize first the vector between the others to have a better projection to the low modes
	/*for (int j = 0; j < multigrid_vector_t::Layout::basisDimension; ++j) {
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
	AlgebraUtils::normalize(vectorspace[index]);*/

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

