#include "TwoFlavorHMCUpdater.h"
#include "dirac_operators/DiracOperator.h"
#include "inverters/BiConjugateGradient.h"
#include "inverters/ConjugateGradient.h"
#include "hmc_integrators/Integrate.h"
#include "algebra_utils/AlgebraUtils.h"
//#define DEBUGFORCE
//#define REVERSIBILITY_CHECK

#ifdef DEBUGFORCE
#include "hmc_forces/TestForce.h"
#endif

namespace Update {

TwoFlavorHMCUpdater::TwoFlavorHMCUpdater() : LatticeSweep(), action(0), gaugeAction(0), fermionAction(0), diracOperator(0), solver(0) { }

TwoFlavorHMCUpdater::TwoFlavorHMCUpdater(const TwoFlavorHMCUpdater& toCopy) : LatticeSweep(toCopy), action(0), gaugeAction(0), fermionAction(0), diracOperator(0), solver(0) { }

TwoFlavorHMCUpdater::~TwoFlavorHMCUpdater() {
	delete action;
	delete solver;
}

void TwoFlavorHMCUpdater::execute(environment_t& environment) {
	//Get the solver for the inversions
	if (solver == 0) solver = new BiConjugateGradient();

	//Initialize the momenta
	this->randomMomenta(momenta);
	//Copy the lattice and the environment
	environmentNew = environment;
	//Initialize the pseudofermion field
	this->generateGaussianDiracVector(tmp_pseudofermion);
	long_real_t oldPseudoFermionEnergy = AlgebraUtils::squaredNorm(tmp_pseudofermion);

	if (diracOperator == 0) diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());

	//Heat bath by multiply the tmp_pseudofermion with the dirac operator
	diracOperator->multiply(pseudofermion, tmp_pseudofermion);

	//Get the gauge action
	if (gaugeAction == 0) gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

	//Set the dirac operator to the new environment
	diracOperator->setLattice(environmentNew.getFermionLattice());
	//Get the fermion action
	if (fermionAction == 0) {
		fermionAction = new TwoFlavorFermionAction(diracOperator);
	}
	else {
		fermionAction->setDiracOperator(diracOperator);
	}

	fermionAction->setPseudoFermion(&pseudofermion);
	fermionAction->setForcePrecision(environment.configurations.get<double>("force_inverter_precision"));

	//Get the global action
	if (action == 0) {
		action = new TwoFlavorAction(gaugeAction,fermionAction);
	}

#ifdef DEBUGFORCE
	//Test of the force
	TestForce testForce;
	extended_gauge_lattice_t tmp;
	//Internal calculations needed
	environmentNew.gaugeLinkConfiguration.updateHalo();
	environmentNew.synchronize();
	action->updateForce(tmp, environmentNew);
	testForce.genericTestForce(environmentNew, action, tmp[5][2], 5, 2);
#endif

	//Get the initial energy of momenta
	long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
	//Get the initial energy of the lattice
	long_real_t oldLatticeEnergy = gaugeAction->energy(environment);

	//The t-length of a single integration step
	real_t t_length = environment.configurations.get<double>("hmc_t_length");
	std::vector<Force*> forces;
	//The numbers of integration steps
	std::vector<unsigned int> numbers_steps = environment.configurations.get< std::vector<unsigned int> >("number_hmc_steps");
	if (numbers_steps.size() == 1) {
		forces.push_back(action);
	} else if (numbers_steps.size() == 2) {
		forces.push_back(fermionAction);
		forces.push_back(gaugeAction);
	}
	else {
		std::cout << "TwoFlavorHMCUpdater::Warning, two flavor does not support more than two multiple time integration!" << std::endl;
		numbers_steps.resize(1);
		forces.push_back(action);
	}
	Integrate* integrate = Integrate::getInstance(environment.configurations.get<std::string>("name_integrator"));

	//Integrate numerically the equation of motion
	integrate->integrate(environmentNew, momenta, forces, numbers_steps, t_length);
#ifdef REVERSIBILITY_CHECK
	integrate->integrate(environmentNew, momenta, forces, numbers_steps, -t_length);
#endif

	//Get the final energy of momenta
	long_real_t newMomentaEnergy = this->momentaEnergy(momenta);
	//Get the final energy of the lattice
	long_real_t newLatticeEnergy = gaugeAction->energy(environmentNew);
	//Get the final energy of the pseudofermions
	diracOperator->setLattice(environmentNew.getFermionLattice());
	solver->setPrecision(environment.configurations.get<double>("metropolis_inverter_precision"));
	solver->setMaximumSteps(environment.configurations.get<unsigned int>("metropolis_inverter_max_steps"));
	solver->solve(diracOperator,pseudofermion,tmp_pseudofermion);
	long_real_t newPseudoFermionEnergy = AlgebraUtils::squaredNorm(tmp_pseudofermion);

	//action->setPseudoFermion(&pseudofermion);TODO why?

	//Global Metropolis Step
	bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy, newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy);
	if (metropolis) {
		environment = environmentNew;
	}

	delete integrate;
}

} /* namespace Update */
