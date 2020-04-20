#include "SquareOverlapOperator.h"
#include "hmc_forces/OverlapFermionForce.h"

namespace Update {

SquareOverlapOperator::SquareOverlapOperator() : DiracOperator(), overlapOperator(0) { }

SquareOverlapOperator::SquareOverlapOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), overlapOperator(0) { }

SquareOverlapOperator::~SquareOverlapOperator() { }

void SquareOverlapOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	overlapOperator->setGamma5(gamma5);
	overlapOperator->multiply(tmp, input);
	overlapOperator->multiply(output, tmp);
}

void SquareOverlapOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	overlapOperator->setGamma5(gamma5);
	overlapOperator->multiply(tmp, vector1);
	overlapOperator->multiplyAdd(output, tmp, vector2, alpha);
}

void SquareOverlapOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	overlapOperator->setKappa(_kappa);
}

void SquareOverlapOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	overlapOperator->setLattice(_lattice);
}

void SquareOverlapOperator::setMass(real_t _mass) {
	overlapOperator->setMass(_mass);
}

real_t SquareOverlapOperator::getMass() const {
	return overlapOperator->getMass();
}

void SquareOverlapOperator::setOverlapOperator(OverlapOperator* _overlapOperator) {
	overlapOperator = _overlapOperator;
}

OverlapOperator* SquareOverlapOperator::getOverlapOperator() const {
	return overlapOperator;
}

void SquareOverlapOperator::setSquareRootApproximation(const Polynomial& _squareRootApproximation) {
	overlapOperator->setSquareRootApproximation(_squareRootApproximation);
}

Polynomial& SquareOverlapOperator::getSquareRootApproximation() {
	return overlapOperator->getSquareRootApproximation();
}

const Polynomial& SquareOverlapOperator::getSquareRootApproximation() const {
	return overlapOperator->getSquareRootApproximation();
}

FermionForce* SquareOverlapOperator::getForce() const {
	std::cout << "Gauge force not implemented for SquareOverlapOperator, return that for OverlapOperator" << std::endl;
	return overlapOperator->getForce();
}

} /* namespace Update */
