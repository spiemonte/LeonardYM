#include "SquareComplementBlockDiracWilsonOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

SquareComplementBlockDiracWilsonOperator::SquareComplementBlockDiracWilsonOperator() : DiracOperator(), diracWilsonOperator(), blockDiracWilsonOperator(), squareBlockDiracWilsonOperator(), log(false), counterSteps(0) { }

SquareComplementBlockDiracWilsonOperator::SquareComplementBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : DiracOperator(_lattice, _kappa), diracWilsonOperator(_lattice, _kappa), blockDiracWilsonOperator(_lattice, _kappa), squareBlockDiracWilsonOperator(_lattice, _kappa), log(false), counterSteps(0) { }

SquareComplementBlockDiracWilsonOperator::~SquareComplementBlockDiracWilsonOperator() { }

void SquareComplementBlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	//First we apply D
	diracWilsonOperator.multiply(output, input);
	
	//Now we apply D_1^(-2) to (D.input)
	biConjugateGradient.solve(&squareBlockDiracWilsonOperator, output, tmpVector);
	counterSteps += biConjugateGradient.getLastSteps();	
	
	//Then we apply again D
	diracWilsonOperator.multiply(output, tmpVector);
	
	if (log && isOutputProcess()) std::cout << "ComplementBlockDiracWilsonOperator::Complement inner inversion done in: " << biConjugateGradient.getLastSteps() << std::endl;
}

void SquareComplementBlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	//First we apply D
	diracWilsonOperator.multiply(output, vector1);
	
	//Now we apply D_1^(-2) to (D.input)
	biConjugateGradient.solve(&squareBlockDiracWilsonOperator, output, tmpVector);
	counterSteps += biConjugateGradient.getLastSteps();	
	
	//Then we apply again D
	diracWilsonOperator.multiply(output, tmpVector);

	//Now we add the shift
	for (int site = 0; site < output.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] += alpha*vector2[site][mu];
	}
	
	if (log && isOutputProcess()) std::cout << "ComplementBlockDiracWilsonOperator::Complement inner inversion done in: " << biConjugateGradient.getLastSteps() << std::endl;
}

void SquareComplementBlockDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	blockDiracWilsonOperator.setLattice(_lattice);
	squareBlockDiracWilsonOperator.setLattice(_lattice);
	diracWilsonOperator.setLattice(_lattice);
}

void SquareComplementBlockDiracWilsonOperator::setBlockSize(const std::vector<unsigned int>& _blockSize) {
	blockDiracWilsonOperator.setBlockSize(_blockSize);
	squareBlockDiracWilsonOperator.setBlockSize(_blockSize);
}

std::vector<unsigned int> SquareComplementBlockDiracWilsonOperator::getBlockSize() const {
	return blockDiracWilsonOperator.getBlockSize();
}

FermionForce* SquareComplementBlockDiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

void SquareComplementBlockDiracWilsonOperator::setPrecision(const real_t& _precision) {
	biConjugateGradient.setPrecision(_precision);
}

real_t SquareComplementBlockDiracWilsonOperator::getPrecision() const {
	return biConjugateGradient.getPrecision();
}

void SquareComplementBlockDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	blockDiracWilsonOperator.setKappa(_kappa);
	diracWilsonOperator.setKappa(_kappa);
	squareBlockDiracWilsonOperator.setKappa(_kappa);
}

} /* namespace Update */
