#include "NFlavorFermionAction.h"
#include "inverters/MultishiftSolver.h"
#include "algebra_utils/AlgebraUtils.h"
#include "hmc_forces/SmearingForce.h"
#include "utils/ToString.h"
namespace Update {

NFlavorFermionAction::NFlavorFermionAction(DiracOperator* _squareDiracOperator, DiracOperator* _diracOperator, const std::vector<RationalApproximation>& _rationalApproximations) : FermionicAction(_diracOperator), squareDiracOperator(_squareDiracOperator), forcePrecision(0.00000000001), maxIterations(5000), rationalApproximations(_rationalApproximations), tmp_pseudofermion(0) {
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
}

NFlavorFermionAction::~NFlavorFermionAction() {
	delete squareDiracOperator;
	delete fermionForce;
}

GaugeGroup NFlavorFermionAction::force(const environment_t& env, int site, int mu) const {
	return fermionForce->force(env, fermionForceLattice[site][mu], site, mu);
}

void NFlavorFermionAction::derivative(const environment_t& env) {
#pragma omp parallel for
	for (int site = 0; site < fermionForceLattice.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(fermionForceLattice[site][mu]);
		}
	}

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
			fermionForce->derivative(fermionForceLattice, env.getFermionLattice(), *x, *y, *weight);
			++x;
			++y;
		}
		++X;
		++Y;
	}
}

void NFlavorFermionAction::updateForce(extended_gauge_lattice_t& forceLattice, const environment_t& env) {
	diracOperator->setLattice(env.getFermionLattice());
	squareDiracOperator->setLattice(env.getFermionLattice());
	fermionForce->setLattice(env.getFermionLattice());
	//Solve the dirac equation for all the pseudofermions
	std::vector< std::vector<extended_dirac_vector_t> >::iterator x = Xs.begin();
	std::vector< std::vector<extended_dirac_vector_t> >::iterator y = Ys.begin();
	std::vector<extended_dirac_vector_t*>::const_iterator pseudofermion = pseudofermions.begin();
	std::vector<RationalApproximation>::iterator i;
	for (i = rationalApproximations.begin(); i != rationalApproximations.end(); ++i) {
		//Solve the dirac equation for all the shifts
		i->getMultishiftSolver()->solve(squareDiracOperator, *(*pseudofermion), *x, i->getBetas());
		std::vector<extended_dirac_vector_t>::const_iterator j;
		std::vector<extended_dirac_vector_t>::iterator k;
		for (j = x->begin(), k = y->begin(); j != x->end(); ++j, ++k) {
			diracOperator->multiply(*k, *j);
		}
		++x;
		++y;
		++pseudofermion;
	}

	//Compute first the link derivative
	this->derivative(env);

	typedef extended_gauge_lattice_t::Layout LT;
	//Set BC for fermions!
	try {
		std::string bc = env.configurations.get<std::string>("boundary_conditions");
		if (bc == "fermion_antiperiodic") {
#pragma omp parallel for
			for (int site = 0; site < fermionForceLattice.completesize; ++site) {
				if (LT::globalIndexT(site) == (LT::glob_t - 1)) fermionForceLattice[site][3] = -fermionForceLattice[site][3];
			}
		}
		else if (bc == "fermion_spatial_antiperiodic") {
#pragma omp parallel for
			for (int site = 0; site < fermionForceLattice.completesize; ++site) {
				if (LT::globalIndexZ(site) == (LT::glob_z - 1)) fermionForceLattice[site][2] = -fermionForceLattice[site][2];
			}
		}
		else if (bc == "fermion_periodic") {
		}
		else if (bc == "open") {
		}
		else {
#pragma omp parallel for
			for (int site = 0; site < fermionForceLattice.completesize; ++site) {
				if (LT::globalIndexT(site) == (LT::glob_t - 1)) fermionForceLattice[site][3] = -fermionForceLattice[site][3];
			}
		}
	}
	catch (NotFoundOption& ex) {
#pragma omp parallel for
		for (int site = 0; site < fermionForceLattice.completesize; ++site) {
			if (LT::globalIndexT(site) == (LT::glob_t - 1)) fermionForceLattice[site][3] = -fermionForceLattice[site][3];
		}
	}
		
	fermionForceLattice.updateHalo();
	
	int smearing_levels = env.configurations.get<unsigned int>("stout_smearing_levels");
	if (smearing_levels == 1) {
		real_t rho = env.configurations.get<real_t>("stout_smearing_rho");
		SmearingForce smearingForce;
		
		smearingForce.force(fermionForceLattice, env.gaugeLinkConfiguration, forceLattice, rho);
	}
	else if (smearing_levels > 1) {
		std::cout << "Smearing levels larger than one not supported!" << std::endl;
		exit(47);
	}
	else {
		//Calculate the force directly
#pragma omp parallel for
		for (int site = 0; site < forceLattice.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				forceLattice[site][mu] = this->force(env, site, mu);
			}
		}
		
		forceLattice.updateHalo();//TODO is needed?
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

int NFlavorFermionAction::getForceMaxIterations() const {
	return maxIterations;
}

void NFlavorFermionAction::setForceMaxIterations(int _maxIterations) {
	maxIterations = _maxIterations;
}

} /* namespace Update */
