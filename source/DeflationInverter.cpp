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

namespace Update {

//Projectors as defined in hep_lat:0706.2298
class Projector : public DiracOperator {
public:
	Projector() : DiracOperator() { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		AlgebraUtils::setToZero(output);
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
			std::complex<real_t> scalar_product = static_cast< std::complex<real_t> >(AlgebraUtils::dot(*vectorspace[i], input));
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] += scalar_product*(*vectorspace[i])[site][mu];
				}
			}
		}
		output.updateHalo();
	}

	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = alpha*vector2[site][mu];
			}
		}
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
			std::complex<real_t> scalar_product = static_cast< std::complex<real_t> >(AlgebraUtils::dot(*vectorspace[i], vector1));
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] += scalar_product*(*vectorspace[i])[site][mu];
				}
			}
		}
		output.updateHalo();
	}

	void addVector(reduced_dirac_vector_t* vector) {
		vectorspace.push_back(vector);
	}

	void clearVectorSpace() {
		vectorspace.clear();
	}

	virtual FermionForce* getForce() const {
		return 0;
	}
private:
	//The vector space used for the projection
	std::vector<reduced_dirac_vector_t*> vectorspace;
};

class IdentityMinusProjector : public DiracOperator {
public:
	IdentityMinusProjector() : DiracOperator() { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		output = input;
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
			std::complex<real_t> scalar_product = static_cast< std::complex<real_t> >(AlgebraUtils::dot(*vectorspace[i], input));
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] -= scalar_product*(*vectorspace[i])[site][mu];
				}
			}
		}
		output.updateHalo();
	}

	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
		std::cout << "Not implemented" << std::endl;
		exit(1);
	}

	void addVector(reduced_dirac_vector_t* vector) {
		vectorspace.push_back(vector);
	}

	void clearVectorSpace() {
		vectorspace.clear();
	}

	virtual FermionForce* getForce() const {
		return 0;
	}
private:
	//The vector space used for the projection
	std::vector<reduced_dirac_vector_t*> vectorspace;
};

class LeftProjector : public DiracOperator {
public:
	LeftProjector() : DiracOperator(), update(true) { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		if (update) this->updateLittleOperator();
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = input[site][mu];
			}
		}
		std::complex<real_t> scalar_products[vectorspace.size()];
		std::complex<real_t> projections[vectorspace.size()];
		for (unsigned int l = 0; l < vectorspace.size(); l += 4) {
			//scalar_products[l] = AlgebraUtils::dot(*vectorspace[l], input);
			long_real_t scalar_product0Re = 0, scalar_product0Im = 0, scalar_product1Re = 0, scalar_product1Im = 0, scalar_product2Re = 0, scalar_product2Im = 0, scalar_product3Re = 0, scalar_product3Im = 0;
#pragma omp parallel for reduction(+:scalar_product0Re,scalar_product0Im,scalar_product1Re,scalar_product1Im,scalar_product2Re,scalar_product2Im,scalar_product3Re,scalar_product3Im)
			for (int site = 0; site < input.localsize; ++site) {
				std::complex<real_t> tmp = vector_dot((*vectorspace[l])[site][0], input[site][0]);
				scalar_product0Re += tmp.real();
				scalar_product0Im += tmp.imag();
				tmp = vector_dot((*vectorspace[l+1])[site][1], input[site][1]);
				scalar_product1Re += tmp.real();
				scalar_product1Im += tmp.imag();
				tmp = vector_dot((*vectorspace[l+2])[site][2], input[site][2]);
				scalar_product2Re += tmp.real();
				scalar_product2Im += tmp.imag();
				tmp = vector_dot((*vectorspace[l+3])[site][3], input[site][3]);
				scalar_product3Re += tmp.real();
				scalar_product3Im += tmp.imag();
			}
			reduceAllSum(scalar_product0Re);
			reduceAllSum(scalar_product0Im);
			reduceAllSum(scalar_product1Re);
			reduceAllSum(scalar_product1Im);
			reduceAllSum(scalar_product2Re);
			reduceAllSum(scalar_product2Im);
			reduceAllSum(scalar_product3Re);
			reduceAllSum(scalar_product3Im);
			scalar_products[l+0] = std::complex<real_t>(scalar_product0Re, scalar_product0Im);
			scalar_products[l+1] = std::complex<real_t>(scalar_product1Re, scalar_product1Im);
			scalar_products[l+2] = std::complex<real_t>(scalar_product2Re, scalar_product2Im);
			scalar_products[l+3] = std::complex<real_t>(scalar_product3Re, scalar_product3Im);
		}
		for (unsigned int k = 0; k < vectorspace.size(); ++k) {
			projections[k] = 0.;
			for (unsigned int l = 0; l < vectorspace.size(); ++l) {
				projections[k] += scalar_products[l]*invertedLittleOperator(k,l);
			}
		}
		for (unsigned int k = 0; k < vectorspace.size(); k += 2) {
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						output[site][mu][c] -= (projections[k]*d_dot_vectorspace[k][site][mu][c] + projections[k+1]*d_dot_vectorspace[k+1][site][mu][c]);
					}
				}
			}
		}
		output.updateHalo();
	}

	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
		if (update) this->updateLittleOperator();
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = vector1[site][mu] + alpha*vector2[site][mu];
			}
		}
		for (unsigned int l = 0; l < vectorspace.size(); ++l) {
			std::complex<real_t> scalar_product = static_cast< std::complex<real_t> >(AlgebraUtils::dot(*vectorspace[l], vector1));
			for (unsigned int k = 0; k < vectorspace.size(); ++k) {
#pragma omp parallel for
				for (int site = 0; site < output.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						output[site][mu] -= scalar_product*invertedLittleOperator(k,l)*d_dot_vectorspace[k][site][mu];
					}
				}
			}
		}
		output.updateHalo();
	}

	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
		update = true;
	}

	void addVector(reduced_dirac_vector_t* vector) {
		vectorspace.push_back(vector);
		update = true;
	}

	void clearVectorSpace() {
		vectorspace.clear();
		update = true;
	}

	matrix_t getInvertedDiracOperator() {
		return invertedLittleOperator;
	}

	virtual FermionForce* getForce() const {
		return 0;
	}
private:
	//The projected DiracOperator
	DiracOperator* dirac;
	//The vector space used for the projection
	std::vector<reduced_dirac_vector_t*> vectorspace;
	//The vector space multiplied by D
	std::vector<reduced_dirac_vector_t> d_dot_vectorspace;
	//flag for internal updates
	bool update;

	matrix_t invertedLittleOperator;

	void updateLittleOperator() {
		if (isOutputProcess()) std::cout << "Updating little operator ..." << std::endl;
		d_dot_vectorspace.resize(vectorspace.size());
		for (unsigned int i = 0; i < vectorspace.size(); ++i) {
			dirac->multiply(d_dot_vectorspace[i],*vectorspace[i]);
		}
		matrix_t littleOperator(vectorspace.size(), vectorspace.size());
		for (unsigned int k = 0; k < vectorspace.size(); ++k) {
			for (unsigned int l = 0; l < vectorspace.size(); ++l) {
				littleOperator(k,l) = AlgebraUtils::dot(*vectorspace[k],d_dot_vectorspace[l]);
			}
		}
		invertedLittleOperator = inverse(littleOperator);
		update = false;
	}
};

class RightProjector : public DiracOperator {
public:
	RightProjector() : DiracOperator(), update(true) { }

	virtual void multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
		if (update) this->updateLittleOperator();
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = input[site][mu];
			}
		}
		dirac->multiply(tmp,input);
		std::complex<real_t>* scalar_products = new std::complex<real_t>[vectorspace.size()];
		for (unsigned int l = 0; l < vectorspace.size(); ++l) scalar_products[l] = AlgebraUtils::dot(*vectorspace[l], tmp);
		for (unsigned int k = 0; k < vectorspace.size(); ++k) {
			std::complex<real_t> projection = 0.;
			for (unsigned int l = 0; l < vectorspace.size(); ++l) {
				projection += scalar_products[l]*invertedLittleOperator(k,l);
			}
#pragma omp parallel for
			for (int site = 0; site < output.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] -= projection*(*vectorspace[k])[site][mu];
				}
			}
		}
		output.updateHalo();
	}

	virtual void multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
		if (update) this->updateLittleOperator();
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = vector1[site][mu] + alpha*vector2[site][mu];
			}
		}
		dirac->multiply(tmp,vector1);
		for (unsigned int l = 0; l < vectorspace.size(); ++l) {//TODO simplify
			std::complex<real_t> scalar_product = static_cast< std::complex<real_t> >(AlgebraUtils::dot(*vectorspace[l], tmp));
			for (unsigned int k = 0; k < vectorspace.size(); ++k) {
#pragma omp parallel for
				for (int site = 0; site < output.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						output[site][mu] -= scalar_product*invertedLittleOperator(k,l)*(*vectorspace[k])[site][mu];
					}
				}
			}
		}
		output.updateHalo();
	}

	void setDiracOperator(DiracOperator* _dirac) {
		dirac = _dirac;
		update = true;
	}

	void addVector(reduced_dirac_vector_t* vector) {
		vectorspace.push_back(vector);
		update = true;
	}

	void clearVectorSpace() {
		vectorspace.clear();
		update = true;
	}

	matrix_t getInvertedDiracOperator() {
		return invertedLittleOperator;
	}

	virtual FermionForce* getForce() const {
		return 0;
	}
private:
	//The projected DiracOperator
	DiracOperator* dirac;
	//The vector space used for the projection
	std::vector<reduced_dirac_vector_t*> vectorspace;
	//flag for internal updates
	bool update;
	reduced_dirac_vector_t tmp;

	matrix_t invertedLittleOperator;

	void updateLittleOperator() {
		matrix_t littleOperator(vectorspace.size(), vectorspace.size());
		for (unsigned int l = 0; l < vectorspace.size(); ++l) {//TODO simplify
			dirac->multiply(tmp,*vectorspace[l]);
			for (unsigned int k = 0; k < vectorspace.size(); ++k) {
				littleOperator(k,l) = AlgebraUtils::dot(*vectorspace[k],tmp);
			}
		}
		invertedLittleOperator = inverse(littleOperator);
		update = false;
	}
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
		return 0;
	}
private:
	//The core operator which has to be deflated
	DiracOperator* dirac;
	//The left projector
	LeftProjector* leftProjector;
	reduced_dirac_vector_t tmp;
};

DeflationInverter::DeflationInverter() : localBasis(0), biConjugateGradient(new BiConjugateGradient()), basisDimension(3) { }

DeflationInverter::DeflationInverter(const DeflationInverter& snd) : localBasis(0), biConjugateGradient(new BiConjugateGradient()), basisDimension(snd.basisDimension) { }

DeflationInverter::~DeflationInverter() {
	//delete biConjugateGradient;
}

bool DeflationInverter::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution) {
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution;
	
	lastStep = 0;

	this->generateBasis(dirac);

	std::cout << "DeflationInverter::Using " << totalNumberOfVectors << " vectors as base for deflation" << std::endl;

	LeftProjector* leftProjector = new LeftProjector();
	for (int i = 0; i < totalNumberOfVectors; ++i) leftProjector->addVector(&localBasis[i]);
	leftProjector->setDiracOperator(dirac);
	RightProjector* rightProjector = new RightProjector();
	for (int i = 0; i < totalNumberOfVectors; ++i) rightProjector->addVector(&localBasis[i]);
	rightProjector->setDiracOperator(dirac);
	DeflatedDirac* deflatedDirac = new DeflatedDirac();
	deflatedDirac->setDiracOperator(dirac);
	deflatedDirac->setLeftProjector(leftProjector);

	reduced_dirac_vector_t tmp1,tmp2,tmp3,tmp4;
	dirac->multiply(tmp1,source);
	leftProjector->multiply(tmp2,tmp1);
	rightProjector->multiply(tmp3,source);
	dirac->multiply(tmp4,tmp3);
	std::cout << "Primo test: " << AlgebraUtils::differenceNorm(tmp4,tmp2) << " " << AlgebraUtils::dot(localBasis[0],localBasis[totalNumberOfVectors-1]) << std::endl;

	reduced_dirac_vector_t projectedSource;
	leftProjector->multiply(projectedSource,source);
	if (isOutputProcess()) std::cout << "Starting speed test ... " << std::endl;
	int numberTests = 100;
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);
	for (int i = 0; i < numberTests; ++i) {
		deflatedDirac->multiply(tmp1,source);
	}
	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	if (isOutputProcess()) std::cout << "Timing for DeflatedDiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
	
	biConjugateGradient->setPrecision(precision);
	extended_dirac_vector_t chie, projectedSourcee = projectedSource;
	biConjugateGradient->solve(deflatedDirac,projectedSourcee,chie);//TODO
	reduced_dirac_vector_t chi = chie;

	/*biConjugateGradient->setPrecision(precision);
	reduced_dirac_vector_t chie, projectedSourcee = projectedSource, projections1, projections2;
	IdentityMinusProjector* projector = new IdentityMinusProjector();
	for (int i = 0; i < totalNumberOfVectors; ++i) projector->addVector(&localBasis[i]);
	projector->multiply(projections1, projectedSource);
	biConjugateGradient->solve(dirac,projections1,chie);//TODO
	projector->multiply(projections2,chie);
	reduced_dirac_vector_t chi = projections2;*/
	
	lastStep += biConjugateGradient->getLastSteps();

	std::cout << "Deflation inverter::Convergence in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;

	reduced_dirac_vector_t chiProjected;
	rightProjector->multiply(chiProjected,chi);

	matrix_t invertedLittleOperator = leftProjector->getInvertedDiracOperator();

#pragma omp parallel for
	for (int site = 0; site < solution.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			solution[site][mu] = chiProjected[site][mu];
		}
	}

	for (int k = 0; k < totalNumberOfVectors; ++k) {
		std::complex<real_t> projectionFactor = 0.;
		for (int l = 0; l < totalNumberOfVectors; ++l) {
			projectionFactor += static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[l],source))*invertedLittleOperator(k,l);
		}
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] += projectionFactor*localBasis[k][site][mu];
			}
		}
	}




	original_solution = solution;

	if (localBasis != 0) delete[] localBasis;
	localBasis = 0;
	return true;

}

void DeflationInverter::setPrecision(double _epsilon) {
	biConjugateGradient->setPrecision(_epsilon);
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

void DeflationInverter::setBlockDivision(int _blockDivision) {
	blockDivision = _blockDivision;
}

int DeflationInverter::getBlockDivision() const {
	return blockDivision;
}

void DeflationInverter::generateBasis(DiracOperator* dirac) {

	typedef reduced_dirac_vector_t::Layout LT;
	
	totalNumberOfVectors = 4*basisDimension;
	if (localBasis != 0) delete[] localBasis;
	localBasis = new reduced_dirac_vector_t[totalNumberOfVectors];
	
	biConjugateGradient->setPrecision(0.00001);
	
	for (int i = 0; i < totalNumberOfVectors; i += 4) {
		extended_dirac_vector_t randomVector, tmp;
		AlgebraUtils::setToZero(randomVector);
		biConjugateGradient->solve(dirac,randomVector,tmp);//TODO
		lastStep += biConjugateGradient->getLastSteps();
		reduced_dirac_vector_t base = tmp;
		for (int site = 0; site < LT::localsize; ++site) {
			localBasis[i][site][0] = tmp[site][0];
			set_to_zero(localBasis[i][site][1]);
			set_to_zero(localBasis[i][site][2]);
			set_to_zero(localBasis[i][site][3]);
			localBasis[i].updateHalo();

			set_to_zero(localBasis[i+1][site][0]);
			localBasis[i+1][site][1] = tmp[site][1];
			set_to_zero(localBasis[i+1][site][2]);
			set_to_zero(localBasis[i+1][site][3]);
			localBasis[i+1].updateHalo();

			set_to_zero(localBasis[i+2][site][0]);
			set_to_zero(localBasis[i+2][site][1]);
			localBasis[i+2][site][2] = tmp[site][2];
			set_to_zero(localBasis[i+2][site][3]);
			localBasis[i+2].updateHalo();

			set_to_zero(localBasis[i+3][site][0]);
			set_to_zero(localBasis[i+3][site][1]);
			set_to_zero(localBasis[i+3][site][2]);
			localBasis[i+3][site][3] = tmp[site][3];
			localBasis[i+3].updateHalo();
		}
	}
		

	for (int i = 0; i < totalNumberOfVectors; ++i) {
		for (int j = 0; j < i; ++j) {
			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(localBasis[j],localBasis[i]));
#pragma omp parallel for
			for (int site = 0; site < localBasis[i].localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					localBasis[i][site][mu] -= proj*localBasis[j][site][mu];
				}
			}
		}
		localBasis[i].updateHalo();
		AlgebraUtils::normalize(localBasis[i]);
	}

	std::cout << "DeflationInverter::Deflation basis generated in " << lastStep << " inversion steps" << std::endl;

}

} /* namespace Update */
