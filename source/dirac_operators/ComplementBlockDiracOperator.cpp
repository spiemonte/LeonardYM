#include "ComplementBlockDiracOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"
#include "algebra_utils/AlgebraUtils.h"

//#define FULL_LOG

namespace Update {

ComplementBlockDiracOperator::ComplementBlockDiracOperator(DiracOperator* _diracOperator, BlockDiracOperator* _redBlockDiracOperator, BlockDiracOperator* _blackBlockDiracOperator) : BlockDiracOperator(), diracOperator(_diracOperator), redBlockDiracOperator(_redBlockDiracOperator), blackBlockDiracOperator(_blackBlockDiracOperator), gmresr(new GMRESR()) { //biConjugateGradient(new BiConjugateGradient()) {
	//biConjugateGradient->setPrecision(0.00000001);
	gmresr->setPrecision(0.00000001);
}

ComplementBlockDiracOperator::ComplementBlockDiracOperator(const ComplementBlockDiracOperator& toCopy) : BlockDiracOperator(toCopy), diracOperator(toCopy.diracOperator), redBlockDiracOperator(toCopy.redBlockDiracOperator), blackBlockDiracOperator(toCopy.blackBlockDiracOperator), gmresr(new GMRESR()) {
	gmresr->setPrecision(0.00000001);
}

ComplementBlockDiracOperator::~ComplementBlockDiracOperator() {
	//delete biConjugateGradient;
	//biConjugateGradient = 0;
	delete gmresr;
	gmresr = 0;
}

void ComplementBlockDiracOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
#ifdef FULL_LOG
	int steps = 0;
#endif
	
	//biConjugateGradient->solve(redBlockDiracOperator, input, tmp1);
	gmresr->solve(redBlockDiracOperator, input, tmp1);
	tmp1.updateHalo();
	redBlockDiracOperator->project(tmp1);

#ifdef FULL_LOG
	steps += biConjugateGradient->getLastSteps();
#endif
	
	diracOperator->multiply(tmp2,tmp1);
#pragma omp parallel for
	for (int site = 0; site < tmp2.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmp2[site][mu] = input[site][mu] - tmp2[site][mu];
		}
	}
	//biConjugateGradient->solve(blackBlockDiracOperator, tmp2, output);
	gmresr->solve(blackBlockDiracOperator, tmp2, output);
	output.updateHalo();
	blackBlockDiracOperator->project(output);
	
#ifdef FULL_LOG
	steps += biConjugateGradient->getLastSteps();
#endif
	
#pragma omp parallel for
	for (int site = 0; site < tmp2.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			output[site][mu] = output[site][mu] + tmp1[site][mu];
		}
	}
	
#ifdef FULL_LOG
	if (isOutputProcess()) std::cout << "ComplementBlockDiracWilsonOperator::Complement Dirac operator evaluated with " << steps << " steps."<< std::endl;
#endif
}

void ComplementBlockDiracOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	//TODO: to be implemented
	/*biConjugateGradient->solve(redBlockDiracOperator, input, tmp1);
	diracOperator->multiplyAdd(tmp2,tmp1);
#pragma omp parallel for
	for (int site = 0; site < tmp2.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			tmp2[site][mu] = input[site][mu] - tmp2[site][mu];
		}
	}
	biConjugateGradient->solve(blackBlockDiracOperator, tmp2, output);
	
#pragma omp parallel for
	for (int site = 0; site < tmp2.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			output[site][mu] = output[site][mu] + tmp1[site][mu];
		}
	}*/
}

void ComplementBlockDiracOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	diracOperator->setLattice(_lattice);
	redBlockDiracOperator->setLattice(_lattice);
	blackBlockDiracOperator->setLattice(_lattice);
}

FermionForce* ComplementBlockDiracOperator::getForce() const {
	return diracOperator->getForce();
}

void ComplementBlockDiracOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracOperator->setKappa(_kappa);
	redBlockDiracOperator->setKappa(_kappa);
	blackBlockDiracOperator->setKappa(_kappa);
}

void ComplementBlockDiracOperator::setPrecision(double precision) {
	gmresr->setPrecision(precision);
}

void ComplementBlockDiracOperator::setMaximumSteps(int steps) {
	gmresr->setMaximumSteps(steps);
}

} /* namespace Update */
