#include "OverlapOperator.h"
#include "hmc_forces/OverlapFermionForce.h"

namespace Update {

OverlapOperator::OverlapOperator() : DiracOperator(), diracWilsonOperator(), squareDiracWilsonOperator() { }

OverlapOperator::OverlapOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), diracWilsonOperator(_lattice, _kappa, _gamma5), squareDiracWilsonOperator(_lattice, _kappa, _gamma5) { }

OverlapOperator::~OverlapOperator() { }

void OverlapOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	squareDiracWilsonOperator.setGamma5(true);
	diracWilsonOperator.setGamma5(true);
	squareRootApproximation.evaluate(&squareDiracWilsonOperator, tmp1, input);
	diracWilsonOperator.multiply(tmp2, tmp1);
	
	if (gamma5) {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*input[site][0] + (0.5-mass/2.)*tmp2[site][0];
			output[site][1] = (0.5+mass/2.)*input[site][1] + (0.5-mass/2.)*tmp2[site][1];
			output[site][2] = -(0.5+mass/2.)*input[site][2] + (0.5-mass/2.)*tmp2[site][2];
			output[site][3] = -(0.5+mass/2.)*input[site][3] + (0.5-mass/2.)*tmp2[site][3];
		}
	}
	else {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*input[site][0] + (0.5-mass/2.)*tmp2[site][0];
			output[site][1] = (0.5+mass/2.)*input[site][1] + (0.5-mass/2.)*tmp2[site][1];
			output[site][2] = (0.5+mass/2.)*input[site][2] - (0.5-mass/2.)*tmp2[site][2];
			output[site][3] = (0.5+mass/2.)*input[site][3] - (0.5-mass/2.)*tmp2[site][3];
		}
	}

	output.updateHalo();
}

void OverlapOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	squareDiracWilsonOperator.setGamma5(true);
	diracWilsonOperator.setGamma5(true);
	squareRootApproximation.evaluate(&squareDiracWilsonOperator, tmp1, vector1);
	diracWilsonOperator.multiply(tmp2, tmp1);

	if (gamma5) {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*vector1[site][0] + (0.5-mass/2.)*tmp2[site][0] + alpha*vector2[site][0];
			output[site][1] = (0.5+mass/2.)*vector1[site][1] + (0.5-mass/2.)*tmp2[site][1] + alpha*vector2[site][1];
			output[site][2] = -(0.5+mass/2.)*vector1[site][2] + (0.5-mass/2.)*tmp2[site][2] + alpha*vector2[site][2];
			output[site][3] = -(0.5+mass/2.)*vector1[site][3] + (0.5-mass/2.)*tmp2[site][3] + alpha*vector2[site][3];
		}
	}
	else {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*vector1[site][0] + (0.5-mass/2.)*tmp2[site][0] + alpha*vector2[site][0];
			output[site][1] = (0.5+mass/2.)*vector1[site][1] + (0.5-mass/2.)*tmp2[site][1] + alpha*vector2[site][0];
			output[site][2] = (0.5+mass/2.)*vector1[site][2] - (0.5-mass/2.)*tmp2[site][2] + alpha*vector2[site][0];
			output[site][3] = (0.5+mass/2.)*vector1[site][3] - (0.5-mass/2.)*tmp2[site][3] + alpha*vector2[site][0];
		}
	}

	output.updateHalo();
}

void OverlapOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracWilsonOperator.setKappa(_kappa);
	squareDiracWilsonOperator.setKappa(_kappa);
}

void OverlapOperator::setMass(real_t _mass) {
	mass = _mass;
}

real_t OverlapOperator::getMass() const {
	return mass;
}

void OverlapOperator::setSquareRootApproximation(const Polynomial& _squareRootApproximation) {
	squareRootApproximation = _squareRootApproximation;
}

Polynomial& OverlapOperator::getSquareRootApproximation() {
	return squareRootApproximation;
}

const Polynomial& OverlapOperator::getSquareRootApproximation() const {
	return squareRootApproximation;
}

void OverlapOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	diracWilsonOperator.setLattice(_lattice);
	squareDiracWilsonOperator.setLattice(_lattice);
}

FermionForce* OverlapOperator::getForce() const {
	return new OverlapFermionForce(kappa, mass, &squareRootApproximation);
}

} /* namespace Update */
