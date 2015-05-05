/*
 * DeflationInverter.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "DeflationInverter.h"
#include "BiConjugateGradient.h"
#include "AlgebraUtils.h"
#include "ToString.h"
#include "ConjugateGradient.h"
//#include "./vectorclass/vectorclass.h"
//#include "./vectorclass/complexvec.h"

const int spin_deflation_size = 1;

namespace Update {

//Projectors as defined in hep_lat:0706.2298
class Projector {
public:
	Projector(reduced_index_lattice_t& _blockIndex, int _totalNumberOfBlocks) : update(true), blockIndex(_blockIndex), totalNumberOfBlocks(_totalNumberOfBlocks) { }

	//To be implemented by the specific projector
	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) = 0;

	//To set the Dirac operator
	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
		update = true;
	}

	//Add a vector to the vector space used for block deflation
	void addVector(const reduced_dirac_vector_t& vector) {
		vectorspace.push_back(vector);
		update = true;
	}

	//Remove all the vectors in the vector space
	void clearVectorSpace() {
		vectorspace.clear();
		update = true;
	}

	//Get the inverse of the little dirac operator
	matrix_t getInverseDiracOperator() {
		return inverseLittleOperator;
	}

	//Set the inverse of the little dirac operator
	void setInverseDiracOperator(matrix_t _inverseLittleOperator) {
		inverseLittleOperator = _inverseLittleOperator;
		update = false;
	}

protected:
	//The projected DiracOperator
	DiracOperator* dirac;
	//The vector space used for the projection
	std::vector<reduced_dirac_vector_t> vectorspace;
	//Temporary vectors to store a single vector in the vector space multiplied by D
	reduced_dirac_vector_t DBlockProjected;
	//Temporary vector to project only a block
	reduced_dirac_vector_t blockProjected;
	//flag for internal updates
	bool update;
	//The lattice with the block indeces
	reduced_index_lattice_t blockIndex;
	//The total number of blocks used
	int totalNumberOfBlocks;

	//The little operator and its inverse
	matrix_t inverseLittleOperator;
	matrix_t littleOperator;

public:

	void updateLittleOperator() {
		if (isOutputProcess()) std::cout << "Projector::Updating little operator ..." << std::endl;
		
		littleOperator.resize(totalNumberOfBlocks*vectorspace.size(), totalNumberOfBlocks*vectorspace.size());
		
		//For the vectors in the basis
		for (unsigned int k1 = 0; k1 < vectorspace.size(); ++k1) {
			//For all the blocks
			for (int i1 = 0; i1 < totalNumberOfBlocks; ++i1) {
				//We extract a single vector which is zero always except in the selected block
#pragma omp parallel for
				for (int site = 0; site < vectorspace[k1].completesize; ++site) {
					if (blockIndex[site] == i1) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							for (int c = 0; c < diracVectorLength; ++c) blockProjected[site][mu][c] = vectorspace[k1][site][mu][c];
						}
					}
					else {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							for (int c = 0; c < diracVectorLength; ++c) blockProjected[site][mu][c] = 0;
						}
					}
				}
			
				//We multiply it by the dirac operator
				dirac->multiply(DBlockProjected,blockProjected);
				
#ifndef MULTITHREADING
				int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
				int numberProcessors = omp_get_max_threads();
#endif
				
				//We compute the projections with all the other elements in the basis
				//TODO we use hermitian condition (if available) to reduce the calculations
				//It does not work if D is not hermitian
				for (unsigned int k2 = k1; k2 < vectorspace.size(); ++k2) {
					//We initialize and set to zero the projections
					std::complex<real_t> proj[totalNumberOfBlocks][numberProcessors];
					for (int k = 0; k < totalNumberOfBlocks; ++k) {
						for (int p = 0; p < numberProcessors; ++p) proj[k][p] = 0.;
					}

					//Compute the scalar product
#pragma omp parallel for
					for (int site = 0; site < vectorspace[k2].localsize; ++site) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
								proj[blockIndex[site]][0] += conj(vectorspace[k2][site][mu][c])*DBlockProjected[site][mu][c];
#endif
#ifdef MULTITHREADING
								proj[blockIndex[site]][omp_get_thread_num()] += conj(vectorspace[k2][site][mu][c])*DBlockProjected[site][mu][c];
#endif
							}
						}
					}
					
					//Collect the results from threads and mpi processes
					for (int k = 0; k < totalNumberOfBlocks; ++k) {
						for (int p = 1; p < numberProcessors; ++p) proj[k][0] += proj[k][p];
					}
					for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[k][0]);						
					
					//Set the little dirac operator
					for (int i2 = 0; i2 < totalNumberOfBlocks; ++i2) {
						littleOperator(totalNumberOfBlocks*k2 + i2,totalNumberOfBlocks*k1 + i1) = proj[i2][0];
						//Here we assume the hermitian condition for D
						littleOperator(totalNumberOfBlocks*k1 + i1,totalNumberOfBlocks*k2 + i2) = conj(proj[i2][0]);
					}
				}
			}
		}
		
		inverseLittleOperator = inverse(littleOperator);
		
		update = false;
	}
};

class InverseLittleOperator : public Projector {
public:
	InverseLittleOperator(reduced_index_lattice_t& _blockIndex, int _totalNumberOfBlocks) : Projector(_blockIndex, _totalNumberOfBlocks) { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		if (update) this->updateLittleOperator();

#ifndef MULTITHREADING
		int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
		int numberProcessors = omp_get_max_threads();
#endif
		//Initialize and set to zero the projections
		std::complex<real_t> proj[vectorspace.size()][totalNumberOfBlocks][numberProcessors];
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[m][k][p] = 0.;
			}
		}

		//Compute the projections
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[i][blockIndex[site]][0] += conj(vectorspace[i][site][mu][c])*input[site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[i][blockIndex[site]][omp_get_thread_num()] += conj(vectorspace[i][site][mu][c])*input[site][mu][c];
#endif
					}
				}
			}
		}

		//Collect the results from threads and mpi processes
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[m][k][0] += proj[m][k][p];
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[m][k][0]);
		}


		//Compute littleDirac*projection
		std::complex<real_t> matrix_factors[vectorspace.size()][totalNumberOfBlocks];
		for (unsigned int m1 = 0; m1 < vectorspace.size(); ++m1) {
			for (int k1 = 0; k1 < totalNumberOfBlocks; ++k1) {
				matrix_factors[m1][k1] = 0.;
				for (unsigned int m2 = 0; m2 < vectorspace.size(); ++m2) {
					for (int k2 = 0; k2 < totalNumberOfBlocks; ++k2) {
						matrix_factors[m1][k1] += proj[m2][k2][0]*inverseLittleOperator(totalNumberOfBlocks*m1 + k1,totalNumberOfBlocks*m2 + k2);
					}
				}
			}
		}
		
		AlgebraUtils::setToZero(output);
		
		//Compute littleDirac*projection*input
		//That's the final result!!!
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						output[site][mu][c] += matrix_factors[m][blockIndex[site]]*vectorspace[m][site][mu][c];
					}
				}
			}
		}
		output.updateHalo();
	}
};


class LittleOperator : public Projector {
public:
	LittleOperator(reduced_index_lattice_t& _blockIndex, int _totalNumberOfBlocks) : Projector(_blockIndex, _totalNumberOfBlocks) { }

	virtual void multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
		if (update) this->updateLittleOperator();

#ifndef MULTITHREADING
		int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
		int numberProcessors = omp_get_max_threads();
#endif
		//Initialize and set to zero the projections
		std::complex<real_t> proj[vectorspace.size()][totalNumberOfBlocks][numberProcessors];
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[m][k][p] = 0.;
			}
		}

		//Compute the projections
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[i][blockIndex[site]][0] += conj(vectorspace[i][site][mu][c])*input[site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[i][blockIndex[site]][omp_get_thread_num()] += conj(vectorspace[i][site][mu][c])*input[site][mu][c];
#endif
					}
				}
			}
		}

		//Collect results from processors and threads
#pragma omp parallel for
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[m][k][0] += proj[m][k][p];
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[m][k][0]);
		}


		//Compute littleDirac*projection
		std::complex<real_t> matrix_factors[vectorspace.size()][totalNumberOfBlocks];
		for (unsigned int m1 = 0; m1 < vectorspace.size(); ++m1) {
			for (int k1 = 0; k1 < totalNumberOfBlocks; ++k1) {
				matrix_factors[m1][k1] = 0.;
				for (unsigned int m2 = 0; m2 < vectorspace.size(); ++m2) {
					for (int k2 = 0; k2 < totalNumberOfBlocks; ++k2) {
						matrix_factors[m1][k1] += proj[m2][k2][0]*littleOperator(totalNumberOfBlocks*m1 + k1,totalNumberOfBlocks*m2 + k2);
					}
				}
			}
		}
		
		AlgebraUtils::setToZero(output);
		
		//Compute littleDirac*projection*input
		//That's the final result!!!
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						output[site][mu][c] += matrix_factors[m][blockIndex[site]]*vectorspace[m][site][mu][c];
					}
				}
			}
		}
		output.updateHalo();
	}
};								

class LeftProjector : public Projector {
public:
	LeftProjector(reduced_index_lattice_t& _blockIndex, int _totalNumberOfBlocks) : Projector(_blockIndex, _totalNumberOfBlocks) { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		if (update) this->updateLittleOperator();

#ifndef MULTITHREADING
		int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
		int numberProcessors = omp_get_max_threads();
#endif
		//Initialize and set to zero the projections
		std::complex<real_t> proj[vectorspace.size()][totalNumberOfBlocks][numberProcessors];
#pragma omp parallel for
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[m][k][p] = 0.;
			}
		}
		
		//Compute the projections
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
#ifndef MULTITHREADING
				int processor = 0;
#endif
#ifdef MULTITHREADING
				int processor = omp_get_thread_num();
#endif		
				int index = blockIndex[site];
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						proj[i][index][processor] += conj(vectorspace[i][site][mu][c])*input[site][mu][c];
					}
				}
			}
		}
		
		/*//Compute the projections
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i]->localsize; ++site) {
#ifndef MULTITHREADING
				int processor = 0;
#endif
#ifdef MULTITHREADING
				int processor = omp_get_thread_num();
#endif		
				int index = blockIndex[site];
				Complex4d a;
				Complex4d b;
				Complex4d result(0.,0.);
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][0][0]));
				b.load(reinterpret_cast<double const*>(&input[site][0][0]));
				result += (~a)*b;
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][0][2]));
				b.load(reinterpret_cast<double const*>(&input[site][0][2]));
				result += (~a)*b;
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][1][1]));
				b.load(reinterpret_cast<double const*>(&input[site][1][1]));
				result += (~a)*b;
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][2][0]));
				b.load(reinterpret_cast<double const*>(&input[site][2][0]));
				result += (~a)*b;
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][2][2]));
				b.load(reinterpret_cast<double const*>(&input[site][2][2]));
				result += (~a)*b;
				a.load(reinterpret_cast<double const*>(&(*vectorspace[i])[site][3][1]));
				b.load(reinterpret_cast<double const*>(&input[site][3][1]));
				result += (~a)*b;
				double final[4];
				result.store(&final[0]);
				proj[i][index][processor] += std::complex<real_t>(final[0]+final[2], final[1]+final[3]);
			}
		}*/

		vector_t projections(vectorspace.size()*totalNumberOfBlocks);
		//Collect the results from threads and mpi processes
#pragma omp parallel for
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				projections[totalNumberOfBlocks*m + k] = 0.;
				for (int p = 0; p < numberProcessors; ++p) projections[totalNumberOfBlocks*m + k] += proj[m][k][p];
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(projections[totalNumberOfBlocks*m + k]);
		}


		//Compute littleDirac*projection
		//std::complex<real_t> matrix_factors[vectorspace.size()][totalNumberOfBlocks];
		vector_t matrix_factors = inverseLittleOperator*projections;
		
		//Set to zero the matrix factors
		/*unsigned int m1;
#pragma omp parallel shared(matrix_factors,proj) private(m1)
		{
#pragma omp for schedule(static)
			for (m1 = 0; m1 < vectorspace.size(); ++m1) {
				for (unsigned int m2 = 0; m2 < vectorspace.size(); ++m2) {
					for (int k1 = 0; k1 < totalNumberOfBlocks; ++k1) {
						for (int k2 = 0; k2 < totalNumberOfBlocks; ++k2) {
							matrix_factors[m1][k1] += proj[m2][k2][0]*inverseLittleOperator(totalNumberOfBlocks*m1 + k1,totalNumberOfBlocks*m2 + k2);
						}
					}
				}
			}
		}*/
		
		AlgebraUtils::setToZero(tmp);

		//Compute littleDirac*projection*input and use a temporary vector
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
#pragma omp parallel for default(shared) schedule(static)
			for (int site = 0; site < output.completesize; ++site) {
				std::complex<real_t> matrix_factor = matrix_factors(totalNumberOfBlocks*m + blockIndex[site]);
				//__builtin_prefetch(&matrix_factors[m][blockIndex[site+1]]);
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						tmp[site][mu][c] -= matrix_factor*vectorspace[m][site][mu][c];
					}
				}
			}
		}
		//output.updateHalo();TODO maybe not needed
		
		
		/*for (unsigned int m = 0; m < vectorspace.size(); ++m) {
#pragma omp parallel for
			for (int site = 0; site < output.completesize; ++site) {
				Complex4d factor(matrix_factors[m][blockIndex[site]].real(),matrix_factors[m][blockIndex[site]].imag());
				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][0][0]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][0][0]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][0][0]));
				}

				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][0][2]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][0][2]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][0][2]));
				}

				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][1][1]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][1][1]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][1][1]));
				}

				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][2][0]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][2][0]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][2][0]));
				}

				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][2][2]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][2][2]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][2][2]));
				}

				{
					Complex4d a;
					a.load(reinterpret_cast<double*>(&(*vectorspace[m])[site][3][1]));
					Complex4d res = a*factor;
					a.load(reinterpret_cast<double*>(&tmp[site][3][1]));
					a -= res;
					a.store(reinterpret_cast<double*>(&tmp[site][3][1]));
				}				
			}
		}*/

		//Multiply the result by the Dirac operator 
		dirac->multiply(output,tmp);

		//The final result!!!
#pragma omp parallel for
		for (int site = 0; site < output.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += input[site][mu];
			}
		}
		//output.updateHalo();TODO maybe not needed
	}

	
private:
	//Temporary vector
	reduced_dirac_vector_t tmp;
};

class RightProjector : public Projector {
public:
	RightProjector(reduced_index_lattice_t& _blockIndex, int _totalNumberOfBlocks) : Projector(_blockIndex, _totalNumberOfBlocks) { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		if (update) this->updateLittleOperator();//TODO do not do it two time!

		//Multiply the result by the Dirac operator and use tmp as temporary vector
		dirac->multiply(tmp,input);		
		
#ifndef MULTITHREADING
		int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
		int numberProcessors = omp_get_max_threads();
#endif
		//Initialize and set to zero the projections
		std::complex<real_t> proj[vectorspace.size()][totalNumberOfBlocks][numberProcessors];
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[m][k][p] = 0.;
			}
		}

		//Compute the projections
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
#pragma omp parallel for
			for (int site = 0; site < vectorspace[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[i][blockIndex[site]][0] += conj(vectorspace[i][site][mu][c])*tmp[site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[i][blockIndex[site]][omp_get_thread_num()] += conj(vectorspace[i][site][mu][c])*tmp[site][mu][c];
#endif
					}
				}
			}
		}

		//Collect the results from threads and mpi processes
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[m][k][0] += proj[m][k][p];
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[m][k][0]);
		}


		//Compute littleDirac*projection
		std::complex<real_t> matrix_factors[vectorspace.size()][totalNumberOfBlocks];
		for (unsigned int m1 = 0; m1 < vectorspace.size(); ++m1) {
			for (int k1 = 0; k1 < totalNumberOfBlocks; ++k1) {
				matrix_factors[m1][k1] = 0.;
				for (unsigned int m2 = 0; m2 < vectorspace.size(); ++m2) {
					for (int k2 = 0; k2 < totalNumberOfBlocks; ++k2) {
						matrix_factors[m1][k1] += proj[m2][k2][0]*inverseLittleOperator(totalNumberOfBlocks*m1 + k1,totalNumberOfBlocks*m2 + k2);
					}
				}
			}
		}
		
		AlgebraUtils::setToZero(output);
		
		//Compute littleDirac*projection*D.input
		for (unsigned int m = 0; m < vectorspace.size(); ++m) {
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						output[site][mu][c] -= matrix_factors[m][blockIndex[site]]*vectorspace[m][site][mu][c];
					}
				}
			}
		}
		output.updateHalo();

		//The final result!!!
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += input[site][mu];
			}
		}
		output.updateHalo();//TODO maybe not needed
	}

private:
	//Temporary vector
	reduced_dirac_vector_t tmp;
};

class DeflatedDirac : public DiracOperator {
public:
	DeflatedDirac() : DiracOperator() { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		dirac->multiply(tmp,input);
		leftProjector->multiply(output,tmp);
	}

	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
		dirac->multiplyAdd(tmp,vector1,vector2,alpha);
		leftProjector->multiply(output,tmp);
	}

	void setLeftProjector(LeftProjector* _leftProjector) {
		leftProjector = _leftProjector;
	}

	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
	}

	virtual FermionForce* getForce() const {
		return dirac->getForce();
	}
private:
	//The core operator which has to be deflated
	DiracOperator* dirac;
	//The left projector
	LeftProjector* leftProjector;
	reduced_dirac_vector_t tmp;
};

DeflationInverter::DeflationInverter(int _recursion) : localBasis(0), conjugateGradient(new ConjugateGradient()), subDeflationInverter(0), leftProjector(0), rightProjector(0), inverseLittleOperator(0), basisDimension(3), recursion(_recursion) { }

DeflationInverter::DeflationInverter(const DeflationInverter& snd) : localBasis(0), conjugateGradient(new ConjugateGradient()), subDeflationInverter(0), leftProjector(0), rightProjector(0), inverseLittleOperator(0), basisDimension(snd.basisDimension),
																		blockSizeX(snd.blockSizeX), blockSizeY(snd.blockSizeY), blockSizeZ(snd.blockSizeZ), blockSizeT(snd.blockSizeT), recursion(snd.recursion) { }

DeflationInverter::~DeflationInverter() {
	if (localBasis != 0) delete[] localBasis;
	if (subDeflationInverter != 0) delete subDeflationInverter;
	localBasis = 0;
	delete conjugateGradient;
}

bool DeflationInverter::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution) {
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution;
	
	lastStep = 0;

	//this->generateBasis(dirac);TODO, the set of the basis is left to the user

	std::cout << "DeflationInverter::Using " << basisDimension << " vectors and " << totalNumberOfBlocks << " blocks as base for deflation" << std::endl;
	
	DeflatedDirac* deflatedDirac = new DeflatedDirac();
	deflatedDirac->setDiracOperator(dirac);
	deflatedDirac->setLeftProjector(leftProjector);
	

	/*reduced_dirac_vector_t tmp1,tmp2,tmp3,tmp4;
	dirac->multiply(tmp1,source);
	leftProjector->multiply(tmp2,tmp1);
	rightProjector->multiply(tmp3,source);
	dirac->multiply(tmp4,tmp3);
	std::cout << "Primo test: " << AlgebraUtils::differenceNorm(tmp4,tmp2) << std::endl;

	leftProjector->multiply(tmp1,source);
	leftProjector->multiply(tmp2,tmp1);
	std::cout << "Secondo test: " << AlgebraUtils::differenceNorm(tmp1,tmp2) << std::endl;

	rightProjector->multiply(tmp1,source);
	rightProjector->multiply(tmp2,tmp1);
	std::cout << "Terzo test: " << AlgebraUtils::differenceNorm(tmp1,tmp2) << std::endl;*/
	
	reduced_dirac_vector_t projectedSource;
	leftProjector->multiply(projectedSource,source);
	
	/*if (isOutputProcess()) std::cout << "Starting speed test ... " << std::endl;
	reduced_dirac_vector_t tmp1;
	int numberTests = 1000;
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);
	for (int i = 0; i < numberTests; ++i) {
		deflatedDirac->multiply(tmp1,source);
	}
	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	if (isOutputProcess()) std::cout << "Timing for DeflatedDiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;*/

	extended_dirac_vector_t chie, projectedSourcee = projectedSource;

	if (recursion != 0) {
		subDeflationInverter->solve(deflatedDirac,projectedSourcee,chie);
	}
	else {
		conjugateGradient->setPrecision(precision);
		conjugateGradient->setMaximumSteps(3000);
		//conjugateGradient->solve(deflatedDirac,projectedSourcee,chie);

		lastStep += conjugateGradient->getLastSteps();
		std::cout << "Deflation inverter::Convergence in " << conjugateGradient->getLastSteps() << " steps." << std::endl;
	}
	
	reduced_dirac_vector_t chi = chie;

	/*biConjugateGradient->setPrecision(precision);
	reduced_dirac_vector_t chie, projectedSourcee = projectedSource, projections1, projections2;
	IdentityMinusProjector* projector = new IdentityMinusProjector();
	for (int i = 0; i < totalNumberOfVectors; ++i) projector->addVector(&localBasis[i]);
	projector->multiply(projections1, projectedSource);
	biConjugateGradient->solve(dirac,projections1,chie);//TODO
	projector->multiply(projections2,chie);
	reduced_dirac_vector_t chi = projections2;*/
	
	

	reduced_dirac_vector_t chiProjected;
	rightProjector->multiply(chiProjected,chi);

	/*InverseLittleOperator* inverseLittleOperator = new InverseLittleOperator(blockIndex,totalNumberOfBlocks);
	inverseLittleOperator->setDiracOperator(dirac);
	for (int i = 0; i < basisDimension; ++i) inverseLittleOperator->addVector(&localBasis[i]);*/
	reduced_dirac_vector_t rho;
	inverseLittleOperator->multiply(rho,source);

	/*LittleOperator* littleOperator = new LittleOperator(blockIndex,totalNumberOfBlocks);
	littleOperator->setDiracOperator(dirac);
	for (int i = 0; i < basisDimension; ++i) littleOperator->addVector(&singlePrecisionLocalBasis[i]);
	reduced_dirac_vector_t rho2;
	littleOperator->multiply(rho2,rho);
	std::cout << "Quarto test: " << AlgebraUtils::differenceNorm(rho2,source) << std::endl;*/

	
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			original_solution[site][mu] = rho[site][mu] + chiProjected[site][mu];
		}
	}

	//solution = solution;
	return true;

}

void DeflationInverter::setPrecision(double _epsilon) {
	conjugateGradient->setPrecision(_epsilon);
	precision = _epsilon;
}

double DeflationInverter::getPrecision() const {
	return precision;
}

unsigned int DeflationInverter::getLastSteps() const {
	return lastStep;
}

void DeflationInverter::setBasisDimension(int _basisDimension) {
	basisDimension = _basisDimension;
}

int DeflationInverter::getBasisDimension() const {
	return basisDimension;
}

void DeflationInverter::setBlockSize(int _blockSize) {
	blockSizeX = _blockSize;
	blockSizeY = _blockSize;
	blockSizeZ = _blockSize;
	blockSizeT = _blockSize;
}

void DeflationInverter::setBlockSizeX(int _blockSizeX) {
	blockSizeX = _blockSizeX;
}

void DeflationInverter::setBlockSizeY(int _blockSizeY) {
	blockSizeY = _blockSizeY;
}

void DeflationInverter::setBlockSizeZ(int _blockSizeZ) {
	blockSizeZ = _blockSizeZ;
}

void DeflationInverter::setBlockSizeT(int _blockSizeT) {
	blockSizeT = _blockSizeT;
}

int DeflationInverter::getBlockVolume() const {
	return blockSizeX*blockSizeY*blockSizeZ*blockSizeT;
}

void DeflationInverter::setRecursion(int _recursion) {
	recursion = _recursion;
	if (subDeflationInverter != 0) subDeflationInverter->setRecursion(recursion - 1);
}

void DeflationInverter::generateBasis(DiracOperator* dirac, LeftProjector* upperLeftProjector) {
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);

	//We reset the counter of inversion steps
	lastStep = 0;
	
	typedef reduced_index_lattice_t::Layout LT;
	//number of blocks in the x,y,z,t direction
	int numberBX = LT::glob_x/blockSizeX + ((LT::glob_x % blockSizeX) != 0 ? 1 : 0);
	int numberBY = LT::glob_y/blockSizeY + ((LT::glob_y % blockSizeY) != 0 ? 1 : 0);
	int numberBZ = LT::glob_z/blockSizeZ + ((LT::glob_z % blockSizeZ) != 0 ? 1 : 0);
	int numberBT = LT::glob_t/blockSizeT + ((LT::glob_t % blockSizeT) != 0 ? 1 : 0);

	totalNumberOfBlocks = numberBX*numberBY*numberBZ*numberBT;
	
	for (int site = 0; site < blockIndex.localsize; ++site) {
		int x = LT::globalIndexX(site)/blockSizeX;
		int y = LT::globalIndexY(site)/blockSizeY;
		int z = LT::globalIndexZ(site)/blockSizeZ;
		int t = LT::globalIndexT(site)/blockSizeT;

		blockIndex[site] = numberBT*(numberBZ*(numberBY*x + y) + z) + t;
	}
	blockIndex.updateHalo();
	
	if (localBasis != 0) delete[] localBasis;
	localBasis = new reduced_dirac_vector_t[basisDimension];

	lastStep = 0;
	
	conjugateGradient->setPrecision(0.0000000001);
	conjugateGradient->setMaximumSteps(300);

	reduced_dirac_vector_t zeroVector;
	reduced_dirac_vector_t randomVector;
	AlgebraUtils::setToZero(zeroVector);
	
	for (int i = 0; i < basisDimension; ++i) {
		//We start with a random vector
		AlgebraUtils::generateRandomVector(randomVector);
		
		/*//We orthogonalize first the vectors between themself to have a better projection to the low modes
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[j],randomVector));
#pragma omp parallel for
			for (int site = 0; site < localBasis[i].completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					randomVector[site][mu] -= proj*localBasis[j][site][mu];
				}
			}
		}*/
		//localBasis[i].updateHalo();TODO maybe not needed
		//AlgebraUtils::normalize(randomVector);

		//We give random vector as guess to solve the omogeneous system
		if (upperLeftProjector == 0) {
			conjugateGradient->setMaximumSteps(500);
			conjugateGradient->solve(dirac,zeroVector,localBasis[i],&randomVector);
			
		}
		else {
			conjugateGradient->setMaximumSteps(250);
			reduced_dirac_vector_t tmpc;
			upperLeftProjector->multiply(tmpc, randomVector);
			//conjugateGradient->solve(dirac,randomVector,localBasis[i]);
			conjugateGradient->solve(dirac,tmpc,localBasis[i]/*,&randomVector*/);
		}
		//biConjugateGradient->solve(dirac,zeroVector,randomVector);
		//biConjugateGradient->solve(dirac,randomVector,localBasis[i]);
		lastStep += conjugateGradient->getLastSteps();
	}

	/*//We orthogonalize first the vectors between themself to have a better projection to the low modes
	for (int i = 0; i < basisDimension; ++i) {
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[j],localBasis[i]));
#pragma omp parallel for
			for (int site = 0; site < localBasis[i].completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					localBasis[i][site][mu] -= proj*localBasis[j][site][mu];
				}
			}
		}
		//localBasis[i].updateHalo();TODO maybe not needed
		AlgebraUtils::normalize(localBasis[i]);
	}*/
	
	//Now we perform a block orthogonalization of the local basis
#ifndef MULTITHREADING
	int numberProcessors = 1;
#endif
#ifdef MULTITHREADING
	int numberProcessors = omp_get_max_threads();
#endif
	
	for (int i = 0; i < basisDimension; ++i) {
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj[totalNumberOfBlocks][numberProcessors];
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 0; p < numberProcessors; ++p) proj[k][p] = 0.;
			}
#pragma omp parallel for
			for (int site = 0; site < localBasis[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
						proj[blockIndex[site]][0] += conj(localBasis[j][site][mu][c])*localBasis[i][site][mu][c];
#endif
#ifdef MULTITHREADING
						proj[blockIndex[site]][omp_get_thread_num()] += conj(localBasis[j][site][mu][c])*localBasis[i][site][mu][c];
#endif
					}
				}
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) {
				for (int p = 1; p < numberProcessors; ++p) proj[k][0] += proj[k][p];
			}
			for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(proj[k][0]);
			//= static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[j],localBasis[i]));
#pragma omp parallel for
			for (int site = 0; site < localBasis[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					localBasis[i][site][mu] -= proj[blockIndex[site]][0]*localBasis[j][site][mu];
				}
			}
		}
		
		//Now we block normalize the vector
		real_t norm[totalNumberOfBlocks][numberProcessors];
		for (int k = 0; k < totalNumberOfBlocks; ++k) {
			for (int p = 0; p < numberProcessors; ++p) norm[k][p] = 0.;
		}
#pragma omp parallel for
		for (int site = 0; site < localBasis[i].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
#ifndef MULTITHREADING
					norm[blockIndex[site]][0] += real(conj(localBasis[i][site][mu][c])*localBasis[i][site][mu][c]);
#endif
#ifdef MULTITHREADING
					norm[blockIndex[site]][omp_get_thread_num()] += real(conj(localBasis[i][site][mu][c])*localBasis[i][site][mu][c]);
#endif
				}
			}
		}
		for (int k = 0; k < totalNumberOfBlocks; ++k) {
			for (int p = 1; p < numberProcessors; ++p) norm[k][0] += norm[k][p];
		}
		for (int k = 0; k < totalNumberOfBlocks; ++k) reduceAllSum(norm[k][0]);
		for (int k = 0; k < totalNumberOfBlocks; ++k) norm[k][0] = sqrt(norm[k][0]);
		//= static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[j],localBasis[i]));
#pragma omp parallel for
		for (int site = 0; site < localBasis[i].localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				localBasis[i][site][mu] = localBasis[i][site][mu]/norm[blockIndex[site]][0];
			}
		}		
		
		localBasis[i].updateHalo();
	}

	/*reduced_dirac_vector_t test;
	//We extract a single vector which is zero always except in the selected block
#pragma omp parallel for
	for (int site = 0; site < localBasis[13].completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			if (blockIndex[site][mu] == 11) test[site][mu] = localBasis[13][site][mu];
			else {
				for (int c = 0; c < diracVectorLength; ++c) test[site][mu][c] = 0;
			}
		}
	}
	for (int i = 0; i < basisDimension; ++i) {
		if (i == 13) std::cout << "Eccoci!" << std::endl;
		std::cout << "Vediamo che succede: " << AlgebraUtils::dot(test,localBasis[i]) << std::endl;
	}*/

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	
	if (isOutputProcess()) std::cout << "DeflationInverter::Deflation basis generated in " << lastStep << " inversion steps" << std::endl;
	if (isOutputProcess()) std::cout << "DeflationInverter::   and in " << elapsed << " s." << std::endl;

	clock_gettime(CLOCK_REALTIME, &start);
	
	leftProjector = new LeftProjector(blockIndex,totalNumberOfBlocks);
	for (int i = 0; i < basisDimension; ++i) leftProjector->addVector(localBasis[i]);
	leftProjector->setDiracOperator(dirac);
	leftProjector->updateLittleOperator();
	
	rightProjector = new RightProjector(blockIndex,totalNumberOfBlocks);
	for (int i = 0; i < basisDimension; ++i) rightProjector->addVector(localBasis[i]);
	rightProjector->setDiracOperator(dirac);
	rightProjector->setInverseDiracOperator(leftProjector->getInverseDiracOperator());
	
	inverseLittleOperator = new InverseLittleOperator(blockIndex,totalNumberOfBlocks);
	for (int i = 0; i < basisDimension; ++i) inverseLittleOperator->addVector(localBasis[i]);
	inverseLittleOperator->setDiracOperator(dirac);
	inverseLittleOperator->setInverseDiracOperator(leftProjector->getInverseDiracOperator());

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	if (isOutputProcess()) std::cout << "DeflationInverter::   little operator updated in " << (elapsed) << " s."<< std::endl;


	if (recursion != 0) {
		DeflatedDirac* deflatedDirac = new DeflatedDirac();
		deflatedDirac->setDiracOperator(dirac);
		deflatedDirac->setLeftProjector(leftProjector);
		
		subDeflationInverter = new DeflationInverter();
		subDeflationInverter->setPrecision(precision);
		subDeflationInverter->setBasisDimension(51);
		subDeflationInverter->setBlockSize(4);
		subDeflationInverter->setRecursion(recursion-1);
		subDeflationInverter->generateBasis(deflatedDirac,leftProjector);

		//delete deflatedDirac;
	}

}

} /* namespace Update */
