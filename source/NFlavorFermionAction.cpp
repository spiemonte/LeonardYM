/*
 * NFlavorFermionAction.cpp
 *
 *  Created on: Aug 8, 2012
 *      Author: spiem_01
 */

#include "NFlavorFermionAction.h"
#include "MultishiftSolver.h"
#include "AlgebraUtils.h"
#include "ChronologicalMultishiftSolver.h"
#include "MMMRMultishiftSolver.h"
#include "MEMultishiftSolver.h"

namespace Update {

std::vector< std::vector<extended_dirac_vector_t> > NFlavorFermionAction::Xs;
std::vector< std::vector<extended_dirac_vector_t> > NFlavorFermionAction::Ys;
extended_dirac_vector_t* NFlavorFermionAction::tmp_pseudofermion = 0;

NFlavorFermionAction::NFlavorFermionAction(DiracOperator* _squareDiracOperator, DiracOperator* _diracOperator, const std::vector<RationalApproximation>& _rationalApproximations) : FermionicAction(_diracOperator), squareDiracOperator(_squareDiracOperator), forcePrecision(0.00000000001), rationalApproximations(_rationalApproximations) {
	fermionForce = diracOperator->getForce();
	//Allocate the memory for all the pseudofermions needed for the calculation of the force ( # of vectors = 2*sum(order(rationalApproximations[i]) )
	if (Xs.size() != rationalApproximations.size() || Ys.size() != rationalApproximations.size()) {
		Xs.resize(rationalApproximations.size());
		Ys.resize(rationalApproximations.size());
		std::vector< std::vector<extended_dirac_vector_t> >::iterator x = Xs.begin();
		std::vector< std::vector<extended_dirac_vector_t> >::iterator y = Ys.begin();
		std::vector<RationalApproximation>::const_iterator i;
		for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
			x->resize(i->getAlphas().size());
			y->resize(i->getAlphas().size());

			//Set the vectors to random, needed for the chronological inverter
			std::vector<extended_dirac_vector_t>::iterator xv = x->begin();
			std::vector<extended_dirac_vector_t>::iterator yv = y->begin();
			while (xv != x->end()) {
				AlgebraUtils::generateRandomVector(*xv);
				AlgebraUtils::generateRandomVector(*yv);
				++xv;
				++yv;
			}

			++x;
			++y;
		}
	}
	if (tmp_pseudofermion == 0) tmp_pseudofermion = new extended_dirac_vector_t;

	multishiftSolver = new MMMRMultishiftSolver();
}

NFlavorFermionAction::~NFlavorFermionAction() {
	delete squareDiracOperator;
	delete fermionForce;
	delete multishiftSolver;
}

GaugeGroup NFlavorFermionAction::force(const environment_t& env, int site, int mu) const {
	GaugeGroup force;
	set_to_zero(force);
	std::vector<RationalApproximation>::const_iterator i;
	std::vector< std::vector<extended_dirac_vector_t> >::const_iterator X = Xs.begin();
	std::vector< std::vector<extended_dirac_vector_t> >::const_iterator Y = Ys.begin();
	for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
		//The vector of the weights (alphas)
		std::vector< real_t > weights = i->getAlphas();
		//the pointer to the single weight
		std::vector< real_t >::const_iterator weight;
		//Take the list of the solutions for the single fermion action
		std::vector<extended_dirac_vector_t>::const_iterator x = X->begin();
		std::vector<extended_dirac_vector_t>::const_iterator y = Y->begin();
		for (weight = weights.begin(); weight != weights.end(); ++weight) {
			//Minus sign on the fermion force!
			force -= (*weight)*fermionForce->force(env, fermionForce->derivative(env.getFermionLattice(), *x, *y, site, mu), site, mu);
			++x;
			++y;
		}
		++X;
		++Y;
	}
	return force;
}

void NFlavorFermionAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	diracOperator->setLattice(env.getFermionLattice());
	squareDiracOperator->setLattice(env.getFermionLattice());
	//Take the multi-shift solver
	multishiftSolver->setPrecision(forcePrecision);
	//Solve the dirac equation for all the pseudofermions
	std::vector< std::vector<extended_dirac_vector_t> >::iterator x = Xs.begin();
	std::vector< std::vector<extended_dirac_vector_t> >::iterator y = Ys.begin();
	std::vector<extended_dirac_vector_t*>::const_iterator pseudofermion = pseudofermions.begin();
	std::vector<RationalApproximation>::const_iterator i;
	for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
		//Solve the dirac equation for all the shifts
		multishiftSolver->solve(squareDiracOperator, *(*pseudofermion), *x, i->getBetas());
		std::vector<extended_dirac_vector_t>::const_iterator j;
		std::vector<extended_dirac_vector_t>::iterator k;
		for (j = x->begin(), k = y->begin(); j != x->end(); ++j, ++k) {
			diracOperator->multiply(*k, *j);
		}
		++x;
		++y;
		++pseudofermion;//TODO
	}

	//Calculate the force
#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			forceLattice[site][mu] = this->force(env, site, mu);
		}
	}

/*#pragma omp parallel for
	for (int site = 0; site < forceLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(forceLattice[site][mu]);
		}
	}

	std::vector< std::vector<extended_dirac_vector_t> >::const_iterator X = Xs.begin();
	std::vector< std::vector<extended_dirac_vector_t> >::const_iterator Y = Ys.begin();
	for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
		//The vector of the weights (alphas)
		std::vector< real_t > weights = i->getAlphas();
		//the pointer to the single weight
		std::vector< real_t >::const_iterator weight;
		//Take the list of the solutions for the single fermion action
		std::vector<extended_dirac_vector_t>::const_iterator x = X->begin();
		std::vector<extended_dirac_vector_t>::const_iterator y = Y->begin();
		for (weight = weights.begin(); weight != weights.end(); ++weight) {
			//Minus sign on the fermion force!
#pragma omp parallel for
			for (int site = 0; site < forceLattice.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					forceLattice[site][mu] -= (*weight)*fermionForce->force(env, fermionForce->derivative(env.getFermionLattice(), *x, *y, site, mu), site, mu);
				}
			}
			++x;
			++y;
		}
		++X;
		++Y;
	}*/

	forceLattice.updateHalo();//TODO is needed?
}

long_real_t NFlavorFermionAction::energy(const environment_t& env) {
	diracOperator->setLattice(env.getFermionLattice());
	squareDiracOperator->setLattice(env.getFermionLattice());
	std::vector<RationalApproximation>::iterator i;
	std::vector<extended_dirac_vector_t*>::const_iterator pseudofermion = pseudofermions.begin();
	long_real_t energy = 0.;
	for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
		i->evaluate(squareDiracOperator, *tmp_pseudofermion, **pseudofermion);
		//Dot with the pseudofermions, no norm!
		energy += real(AlgebraUtils::dot(**pseudofermion, *tmp_pseudofermion));
		++pseudofermion;
	}
	//Plus on the pseudofermion energy!
	return energy;
}

void NFlavorFermionAction::addPseudoFermion(extended_dirac_vector_t* _pseudofermion) {
	pseudofermions.push_back(_pseudofermion);
}

void NFlavorFermionAction::cleanPseudoFermions() {
	pseudofermions.clear();
}

std::vector<extended_dirac_vector_t*> NFlavorFermionAction::getPseudoFermion() const {
	return pseudofermions;
}

double NFlavorFermionAction::getForcePrecision() const {
	return forcePrecision;
}

void NFlavorFermionAction::setForcePrecision(double precision) {
	forcePrecision = precision;
}

} /* namespace Update */
