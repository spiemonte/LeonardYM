/*
 * NFlavorQCDUpdater.cpp
 *
 *  Created on: May 2, 2012
 *      Author: spiem_01
 */

#include "NFlavorBlockUpdater.h"
#include "NFlavorQCDAction.h"
#include "AlgebraUtils.h"
#include "MultishiftSolver.h"
#include "RationalApproximation.h"
#include "BiConjugateGradient.h"
#include "DiracEigenSolver.h"
#include "GaugeAction.h"
#include "Integrate.h"
#include "ToString.h"
#include "GlobalOutput.h"
#include "./dirac_operators/BlockDiracOperator.h"
#include "./dirac_operators/SquareComplementBlockDiracWilsonOperator.h"
#include "./dirac_operators/SquareDiracWilsonOperator.h"
#include <iomanip>

//#define DEBUGFORCE
//#define REVERSIBILITY_CHECK

#ifdef DEBUGFORCE
#include "TestForce.h"
#endif

namespace Update {

NFlavorBlockUpdater::NFlavorBlockUpdater() : MultiStepNFlavorQCDUpdater(), biConjugateGradient(0) { }

NFlavorBlockUpdater::NFlavorBlockUpdater(const NFlavorBlockUpdater& toCopy) : MultiStepNFlavorQCDUpdater(toCopy), biConjugateGradient(0) { }

NFlavorBlockUpdater::~NFlavorBlockUpdater() {
	if (biConjugateGradient != 0) delete biConjugateGradient;
	biConjugateGradient = 0;
}

void NFlavorBlockUpdater::execute(environment_t& environment) {	
	//First we initialize the approximations
	this->initializeApproximations(environment);

	//We check the theory that is simulated
	if (environment.iteration == 0 && environment.sweep == 0) {
		this->checkTheory(environment);
	}

	//Initialize the momenta
	this->randomMomenta(momenta);
	//Copy the environment
	environmentNew = environment;

	//Take the Block Dirac Operator
	if (diracOperator == 0)  diracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());

	//Take the Block Square Dirac Operator
	if (squareDiracOperator == 0) squareDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	squareDiracOperator->setLattice(environment.getFermionLattice());

	//Take the gauge action
	if (gaugeAction == 0) gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

	//Initialize the pseudofermion fields
	std::vector<extended_dirac_vector_t>::iterator i;
	std::vector<RationalApproximation>::iterator j = rationalApproximationsHeatBath.begin();
	long_real_t oldPseudoFermionEnergy = 0.;
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		this->generateGaussianDiracVector(tmp_pseudofermion);
		oldPseudoFermionEnergy += AlgebraUtils::squaredNorm(tmp_pseudofermion);
		//Heat bath by multiply the tmp_pseudofermion to the polynomial of the dirac operator
		j->evaluate(squareDiracOperator, *i, tmp_pseudofermion);//TODO
		++j;
	}

	try {
		std::string flag = environment.configurations.get<std::string>("check_rational_approximations");

		if (environment.sweep == 0 && environment.iteration == 0 && environment.measurement && flag == "true") {
			//Now we test the correctness of rational/heatbath
			long_real_t test = 0.;
			//Now we use the better approximation for the metropolis step
			std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
			for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
				//Now we evaluate it with the rational approximation of the inverse
				rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
				test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
				++rational;
			}
			if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Consistency check for the metropolis: " << test - oldPseudoFermionEnergy << std::endl;
			//Now we use the approximation for the force step
			rational = rationalApproximationsForce[0].begin();
			test = 0.;
			for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
				//Now we evaluate it with the rational approximation of the inverse
				rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
				test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
				++rational;
			}
			if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Consistency check for the first level of the force: " << test - oldPseudoFermionEnergy << std::endl;
		}

	} catch (NotFoundOption& ex) {
		if (isOutputProcess() && environment.sweep == 0 && environment.iteration == 0 && environment.measurement) std::cout << "NFlavorQCDUpdater::No consistency check of metropolis/force approximations!" << std::endl;
	}

	//Get the initial energy of momenta
	long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
	//Get the initial energy of the lattice
	long_real_t oldLatticeEnergy = gaugeAction->energy(environment);

	//Set the dirac operator to the new environment
	squareDiracOperator->setLattice(environmentNew.getFermionLattice());
	diracOperator->setLattice(environmentNew.getFermionLattice());

	//Take the fermion action
	if (fermionAction == 0) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		fermionAction = new NFlavorFermionAction*[numberLevels];
		for (int j = 0; j < numberLevels; ++j) {
			fermionAction[j] = new NFlavorFermionAction(squareDiracOperator, diracOperator, rationalApproximationsForce[j]);
			//TODO dirac operator correctly set?
			for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
				//Add the pseudofermions
				fermionAction[j]->addPseudoFermion(&(*i));
			}
			//Set the precision of the inverter
			fermionAction[j]->setForcePrecision(environment.configurations.get<double>("force_inverter_precision"));
			fermionAction[j]->setForceMaxIterations(environment.configurations.get<unsigned int>("force_inverter_max_steps"));
		}
	}

	//Take the global action
	if (nFlavorQCDAction == 0) nFlavorQCDAction = new NFlavorAction(gaugeAction, fermionAction[0]);//TODO, we skip the other forces

	//The t-length of a single integration step
	real_t t_length = environment.configurations.get<double>("hmc_t_length");
	std::vector<Force*> forces;
	//The numbers of integration steps
	std::vector<unsigned int> numbers_steps = environment.configurations.get< std::vector<unsigned int> >("number_hmc_steps");
	if (numbers_steps.size() == 1) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 1) std::cout << "NFlavorBlockUpdater::Warning, with only one time integration only the first level of the force is used!" << std::endl;
		forces.push_back(nFlavorQCDAction);
	} else if (numbers_steps.size() == 2) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 1) std::cout << "NFlavorBlockUpdater::Warning, with only two time integration only the first level of the force is used!" << std::endl;
		forces.push_back(fermionAction[0]);
		forces.push_back(gaugeAction);
	} else if (numbers_steps.size() == 3) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 2) std::cout << "NFlavorBlockUpdater::Warning, with only three time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[1]);
		forces.push_back(fermionAction[0]);
		forces.push_back(gaugeAction);
	}
	else {
		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Warning, NFlavor does not support more than three time integrations!" << std::endl;
		numbers_steps.resize(1);
		forces.push_back(nFlavorQCDAction);
	}

#ifdef DEBUGFORCE
	//Test of the force
	TestForce testForce;
	extended_gauge_lattice_t tmp;
	//Internal calculations needed
	nFlavorQCDAction->updateForce(tmp, environment);
	testForce.genericTestForce(environment, nFlavorQCDAction, tmp[5][2], 5, 2);
#endif

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
	long_real_t newPseudoFermionEnergy = 0.;
	diracOperator->setLattice(environmentNew.getFermionLattice());//TODO TODO TODO
	squareDiracOperator->setLattice(environmentNew.getFermionLattice());
	//Now we use the better approximation for the metropolis step
	std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		//Now we evaluate it with the rational approximation of the inverse
		rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
		newPseudoFermionEnergy += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
		++rational;
	}

	//Global Metropolis Step
	bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy, newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy);
	if (metropolis) {
		environment = environmentNew;

		real_t alpha = environment.configurations.get<real_t>("theory_power_factor");

		//Now we computed the logarithm of an eventual correction factor correction factor
		log_determinant = this->logDeterminant(environmentNew,alpha);

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");

			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy));
			output->write("hmc_history", 1);

			output->pop("hmc_history");

			output->push("log_determinant");

			output->write("log_determinant", - log_determinant);

			output->pop("log_determinant");
		}
	}
	else {
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");

			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy));
			output->write("hmc_history", 0);

			output->pop("hmc_history");

			output->push("log_determinant");

			output->write("log_determinant", - log_determinant);

			output->pop("log_determinant");
		}
	}

	delete integrate;
}

long_real_t NFlavorBlockUpdater::logDeterminant(const environment_t& environment_new, const environment_t& environment_old, real_t alpha) {
	//Now we get the biconjugate grandient solver
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(0.01);//environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment_new.configurations.get<unsigned int>("generic_inverter_max_steps"));

	//Here we store the partial estimates from the noise technique
	std::vector<long_real_t> stochastic_estimates;

	//Now we get the block dirac wilson operator
	DiracOperator* squareBlockDiracNew = BlockDiracOperator::getInstance(environment_new.configurations.get<std::string>("dirac_operator"), 2, environment_new.configurations);
	squareBlockDiracNew->setLattice(environment_new.getFermionLattice());

	DiracOperator* squareDiracNew = DiracOperator::getInstance(environment_new.configurations.get<std::string>("dirac_operator"), 2, environment_new.configurations);
	squareDiracNew->setLattice(environment_new.getFermionLattice());

	//Now we get the block dirac wilson operator
	DiracOperator* squareBlockDiracOld = BlockDiracOperator::getInstance(environment_old.configurations.get<std::string>("dirac_operator"), 2, environment_old.configurations);
	squareBlockDiracOld->setLattice(environment_old.getFermionLattice());

	DiracOperator* squareDiracOld = DiracOperator::getInstance(environment_old.configurations.get<std::string>("dirac_operator"), 2, environment_old.configurations);
	squareDiracOld->setLattice(environment_old.getFermionLattice());

	int steps = environment_new.configurations.get<unsigned int>("number_block_correction_noise_vectors");

	for (int i = 0; i < steps; ++i) {
		AlgebraUtils::generateRandomComplexGaussianVector(temp1);
		biConjugateGradient->solve(squareDiracNew, temp1, temp2);
		squareBlockDiracNew->multiply(result, temp2);

		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant new: " << AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(temp1, result) << std::endl;

		std::complex<long_real_t> new_log_det = AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(temp1, result);

		biConjugateGradient->solve(squareDiracOld, temp1, temp2);
		squareBlockDiracOld->multiply(result, temp2);

		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant old: " << AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(temp1, result) << std::endl;

		std::complex<long_real_t> old_log_det = AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(temp1, result);
		stochastic_estimates.push_back(real(new_log_det - old_log_det));
	}

	long_real_t average = 0., max = stochastic_estimates.front();
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		average += *i;
		if (*i > max) max = *i;
	}
	average = average/stochastic_estimates.size();

	long_real_t sd = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		sd += (*i - average)*(*i - average);
	}
	sd = sqrt(sd/stochastic_estimates.size());

	long_real_t residuum = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		residuum += exp( *i - max );
	}

	delete squareBlockDiracNew;
	delete squareDiracNew;
	delete squareBlockDiracOld;
	delete squareDiracOld;

	if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Difference fermion action " << average << " with standard deviation " << sd << std::endl;
	return (max - log(static_cast<long_real_t>(stochastic_estimates.size())) + log1p(residuum-1))*alpha;
}

long_real_t NFlavorBlockUpdater::logDeterminant(const environment_t& environment, real_t alpha) {
	//Now we get the biconjugate grandient solver
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));

	//Here we store the partial estimates from the noise technique
	std::vector<long_real_t> stochastic_estimates;

	//Now we get the block dirac wilson operator
	DiracOperator* squareBlockDirac = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	squareBlockDirac->setLattice(environment.getFermionLattice());

	DiracOperator* blockDirac = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	blockDirac->setLattice(environment.getFermionLattice());

	DiracOperator* dirac = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	dirac->setLattice(environment.getFermionLattice());

	int steps = environment.configurations.get<unsigned int>("number_block_correction_noise_vectors");

	for (int i = 0; i < steps; ++i) {
		AlgebraUtils::generateRandomComplexGaussianVector(temp1);
		biConjugateGradient->solve(squareBlockDirac, temp1, temp2);
		blockDirac->multiply(temp3, temp2);
		dirac->multiply(result,temp3);

		long_real_t saturno = real(AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(result, result));

		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant: " << saturno << std::endl;


		stochastic_estimates.push_back(saturno);
	}

	long_real_t average = 0., max = stochastic_estimates.front();
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		average += *i;
		if (*i > max) max = *i;
	}

	average = average/stochastic_estimates.size();

	long_real_t sd = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		sd += (*i - average)*(*i - average);
	}
	sd = sqrt(sd/stochastic_estimates.size());

	long_real_t residuum = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		residuum += exp( *i - max );
	}
	//residuum = residuum/;

	delete squareBlockDirac;
	delete dirac;
	delete blockDirac;

	if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant average " << average << " with standard deviation " << sd << std::endl;
	return (max - log(static_cast<long_real_t>(stochastic_estimates.size())) + log1p(residuum-1))*alpha;
}

long_real_t NFlavorBlockUpdater::logDeterminant_test(const environment_t& environment, real_t alpha) {
	//Now we get the biconjugate grandient solver
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(0.000001);//environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));

	//Here we store the partial estimates from the noise technique
	std::vector<long_real_t> stochastic_estimates;

	//Now we get the block dirac wilson operator
	DiracOperator* squareDiracMass = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	squareDiracMass->setKappa(0.19);
	squareDiracMass->setLattice(environment.getFermionLattice());

	DiracOperator* diracMass = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracMass->setKappa(0.19);
	diracMass->setLattice(environment.getFermionLattice());

	DiracOperator* dirac = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	dirac->setLattice(environment.getFermionLattice());

	int steps = environment.configurations.get<unsigned int>("number_twisted_correction_noise_vectors");

	for (int i = 0; i < steps; ++i) {
		AlgebraUtils::generateRandomComplexGaussianVector(temp1);
		biConjugateGradient->solve(squareDiracMass, temp1, temp2);
		diracMass->multiply(temp3, temp2);
		dirac->multiply(result,temp3);

		long_real_t saturno = real(AlgebraUtils::dot(temp1,temp1) - AlgebraUtils::dot(result, result));

		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant: " << saturno << std::endl;


		stochastic_estimates.push_back(saturno);
	}

	long_real_t average = 0., max = stochastic_estimates.front();
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		average += *i;
		if (*i > max) max = *i;
	}

	average = average/stochastic_estimates.size();

	long_real_t sd = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		sd += (*i - average)*(*i - average);
	}
	sd = sqrt(sd/stochastic_estimates.size());

	long_real_t residuum = 0.;
	for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
		residuum += exp( *i - max );
	}
	//residuum = residuum/;

	delete squareDiracMass;
	delete dirac;
	delete diracMass;

	if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Log determinant average " << average << " with standard deviation " << sd << std::endl;
	return (max - log(static_cast<long_real_t>(stochastic_estimates.size())) + log1p(residuum-1))*alpha;
}

} /* namespace Update */
