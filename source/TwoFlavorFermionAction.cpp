/*
 * TwoFlavorFermionAction.cpp
 *
 *  Created on: Aug 8, 2012
 *      Author: spiem_01
 */

#include "TwoFlavorFermionAction.h"
#include "BiConjugateGradient.h"
#include "ConjugateGradient.h"
#include "AlgebraUtils.h"

namespace Update {

TwoFlavorFermionAction::TwoFlavorFermionAction(DiracOperator* _diracOperator) : FermionicAction(_diracOperator) {
	fermionForce = diracOperator->getForce();
}

TwoFlavorFermionAction::~TwoFlavorFermionAction() {
	delete fermionForce;
}

GaugeGroup TwoFlavorFermionAction::force(const environment_t& env, int site, int mu) const {
	//Minus sign on the fermion force!
	return - fermionForce->force(env, fermionForce->derivative(env.getFermionLattice(), X, Y, site, mu), site, mu);
}

void TwoFlavorFermionAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	diracOperator->setLattice(env.getFermionLattice());
	BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();//TODO TODO TODO
	biConjugateGradient->setPrecision(forcePrecision);
	biConjugateGradient->solve(diracOperator,*pseudofermion,Y);
	biConjugateGradient->solve(diracOperator,Y,X);

	//Calculate the force
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			forceLattice[site][mu] = this->force(env, site, mu);
		}
	}

	forceLattice.updateHalo();//TODO is needed?

	delete biConjugateGradient;
}

long_real_t TwoFlavorFermionAction::energy(const environment_t& env) {
	BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
	diracOperator->setLattice(env.getFermionLattice());
	biConjugateGradient->setPrecision(0.0000000000000000001);
	biConjugateGradient->solve(diracOperator,*pseudofermion,Y);
	delete biConjugateGradient;
	//Plus on the pseudofermion energy!
	return +AlgebraUtils::squaredNorm(Y);
}

void TwoFlavorFermionAction::setPseudoFermion(extended_dirac_vector_t* _pseudofermion) {
	pseudofermion = _pseudofermion;
}

extended_dirac_vector_t* TwoFlavorFermionAction::getPseudoFermion() const {
	return pseudofermion;
}

double TwoFlavorFermionAction::getForcePrecision() const {
	return forcePrecision;
}

void TwoFlavorFermionAction::setForcePrecision(double precision) {
	forcePrecision = precision;
}

} /* namespace Update */
