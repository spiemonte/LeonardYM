#include "MultiStepNFlavorUpdater.h"
#include "actions/NFlavorQCDAction.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/MultishiftSolver.h"
#include "dirac_functions/RationalApproximation.h"
#include "inverters/BiConjugateGradient.h"
#include "inverters/MultiGridMEMultishiftSolver.h"
#include "dirac_operators/ComplementBlockDiracOperator.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "actions/GaugeAction.h"
#include "hmc_integrators/Integrate.h"
#include "utils/ToString.h"
#include "io/GlobalOutput.h"
#include "dirac_operators/SquareTwistedDiracOperator.h"
#include <iomanip>

//#define DEBUGFORCE
//#define REVERSIBILITY_CHECK

#ifdef DEBUGFORCE
#include "hmc_forces/TestForce.h"
#endif

namespace Update {

MultiStepNFlavorUpdater::MultiStepNFlavorUpdater() : LatticeSweep(), nFlavorAction(0), gaugeAction(0), fermionAction(0), squareDiracOperatorMetropolis(0), diracOperatorMetropolis(0), squareDiracOperatorForce(0), diracOperatorForce(0), multishiftSolver(0), blackBlockDiracOperator(0), redBlockDiracOperator(0) { }

MultiStepNFlavorUpdater::MultiStepNFlavorUpdater(const MultiStepNFlavorUpdater& toCopy) : LatticeSweep(toCopy), nFlavorAction(0), gaugeAction(0), fermionAction(0), squareDiracOperatorMetropolis(0), diracOperatorMetropolis(0), squareDiracOperatorForce(0), diracOperatorForce(0), multishiftSolver(0) { }

MultiStepNFlavorUpdater::~MultiStepNFlavorUpdater() {
	if (nFlavorAction != 0) delete nFlavorAction;
	if (multishiftSolver != 0) delete multishiftSolver;
	if (blackBlockDiracOperator != 0) delete blackBlockDiracOperator;
	if (redBlockDiracOperator != 0) delete redBlockDiracOperator;
}

void MultiStepNFlavorUpdater::initializeApproximations(environment_t& environment) {
	double twist = environment.configurations.get<double>("MultiStepNFlavorUpdater::twist");
	if (isOutputProcess()) std::cout << "MultiStepNFlavorUpdater::Using twist " << twist << std::endl;	

	if (multishiftSolver == 0) {
		if (environment.configurations.get<std::string>("MultiStepNFlavorUpdater::multigrid") == "true") {
			DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
			diracOperator->setLattice(environment.getFermionLattice());
			diracOperator->setGamma5(false);

			blackBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Black);
			blackBlockDiracOperator->setLattice(environment.getFermionLattice());
			blackBlockDiracOperator->setGamma5(false);
			blackBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("MultiStepNFlavorUpdater::sap_block_size"));
		
			redBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Red);
			redBlockDiracOperator->setLattice(environment.getFermionLattice());
			redBlockDiracOperator->setGamma5(false);
			redBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("MultiStepNFlavorUpdater::sap_block_size"));

			MultiGridMEMultishiftSolver* multishiftMultiGridSolver = new MultiGridMEMultishiftSolver(environment.configurations.get< unsigned int >("MultiStepNFlavorUpdater::multigrid_basis_dimension"), environment.configurations.get< std::vector<unsigned int> >("MultiStepNFlavorUpdater::multigrid_block_size"), blackBlockDiracOperator, redBlockDiracOperator);
			multishiftMultiGridSolver->setSAPIterations(environment.configurations.get<unsigned int>("MultiStepNFlavorUpdater::sap_iterations"));
			multishiftMultiGridSolver->setSAPMaxSteps(environment.configurations.get<unsigned int>("MultiStepNFlavorUpdater::sap_inverter_max_steps"));
			multishiftMultiGridSolver->setSAPPrecision(environment.configurations.get<real_t>("MultiStepNFlavorUpdater::sap_inverter_precision"));
			multishiftMultiGridSolver->setGMRESIterations(environment.configurations.get<unsigned int>("MultiStepNFlavorUpdater::gmres_inverter_max_steps"));
			multishiftMultiGridSolver->setGMRESPrecision(environment.configurations.get<real_t>("MultiStepNFlavorUpdater::gmres_inverter_precision"));

			multishiftMultiGridSolver->initializeBasis(diracOperator);			
			multishiftSolver = multishiftMultiGridSolver;
			
			delete diracOperator;
			if (isOutputProcess()) std::cout << "MultiStepNFlavorUpdater::Using multigrid inverter and SAP preconditioning ..." << std::endl;
		}
		else {
			multishiftSolver = MultishiftSolver::getInstance("minimal_residual");
		}
	}
	else {
		if (environment.configurations.get<std::string>("MultiStepNFlavorUpdater::multigrid") == "true") {
			MultiGridMEMultishiftSolver* multishiftMultiGridSolver = dynamic_cast<MultiGridMEMultishiftSolver*>(multishiftSolver);
			multishiftMultiGridSolver->initializeBasis(diracOperatorMetropolis);
		}
	}
	//First take the rational function approximation for the heatbath step
	if (rationalApproximationsHeatBath.empty()) {
		int numberPseudofermions = environment.configurations.get< unsigned int >("number_pseudofermions");
		for (int i = 1; i <= numberPseudofermions; ++i) {
			std::vector<real_t> rat = environment.configurations.get< std::vector<real_t> >(std::string("heatbath_rational_fraction_")+toString(i));
			RationalApproximation rational(multishiftSolver);
			rational.setAlphas(std::vector<real_t>(rat.begin(), rat.begin() + rat.size()/2));
			rational.setBetas(std::vector<real_t>(rat.begin() + rat.size()/2, rat.end()));
			//We apply the twist
			for (unsigned int k = 0; k < rational.getBetas().size(); ++k) {
				rational.getBetas()[k] += twist;
			}
			rational.setPrecision(environment.configurations.get<double>("metropolis_inverter_precision"));
			rational.setMaximumRecursion(environment.configurations.get<unsigned int>("metropolis_inverter_max_steps"));
			rationalApproximationsHeatBath.push_back(rational);
		}
		pseudofermions.resize(numberPseudofermions);
	}

	//Then take the rational function approximation for the metropolis step
	if (rationalApproximationsMetropolis.empty()) {
		int numberPseudofermions = environment.configurations.get< unsigned int >("number_pseudofermions");
		for (int i = 1; i <= numberPseudofermions; ++i) {
			std::vector<real_t> rat = environment.configurations.get< std::vector<real_t> >(std::string("metropolis_rational_fraction_")+toString(i));
			RationalApproximation rational(multishiftSolver);
			rational.setAlphas(std::vector<real_t>(rat.begin(), rat.begin() + rat.size()/2));
			rational.setBetas(std::vector<real_t>(rat.begin() + rat.size()/2, rat.end()));
			//We apply the twist
			for (unsigned int k = 0; k < rational.getBetas().size(); ++k) {
				rational.getBetas()[k] += twist;
			}
			rational.setPrecision(environment.configurations.get<double>("metropolis_inverter_precision"));
			rational.setMaximumRecursion(environment.configurations.get<unsigned int>("metropolis_inverter_max_steps"));
			rationalApproximationsMetropolis.push_back(rational);
		}
	}

	//We allow the usage of different precision for different level
	int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
	real_t level_precisions[numberLevels];
	for (int i = 1; i <= numberLevels; ++i) {
		level_precisions[i-1] = environment.configurations.get<double>(std::string("force_inverter_precision_level_")+toString(i));
	}

	//Then take the rational function approximation for the force step
	if (rationalApproximationsForce.empty()) {

		for (int i = 1; i <= numberLevels; ++i) {
			int numberPseudofermions = environment.configurations.get< unsigned int >("number_pseudofermions");
			std::vector<RationalApproximation> levelRationaApproximationForce;
			for (int j = 1; j <= numberPseudofermions; ++j) {
				std::vector<real_t> rat = environment.configurations.get< std::vector<real_t> >(std::string("force_rational_fraction_")+toString(j)+"_level_"+toString(i));
				RationalApproximation rational(multishiftSolver);
				rational.setAlphas(std::vector<real_t>(rat.begin(), rat.begin() + rat.size()/2));
				rational.setBetas(std::vector<real_t>(rat.begin() + rat.size()/2, rat.end()));
				//We apply the twist
				for (unsigned int k = 0; k < rational.getBetas().size(); ++k) {
					rational.getBetas()[k] += twist;
				}
				rational.setPrecision(level_precisions[i - 1]);
				rational.setMaximumRecursion(environment.configurations.get<unsigned int>("force_inverter_max_steps"));
				levelRationaApproximationForce.push_back(rational);
			}
			rationalApproximationsForce.push_back(levelRationaApproximationForce);
		}
	}
}

void MultiStepNFlavorUpdater::checkTheory(const environment_t& environment) const {
	double testerForce = 1., testerMetropolis = 1., testerHeatBath = 1.;
	double epsilon = 0.0001;
	int numberPseudofermions = environment.configurations.get< unsigned int >("number_pseudofermions");
	for (int i = 0; i < numberPseudofermions; ++i) {
		testerForce *= rationalApproximationsForce[0][i].evaluate(2.).real();
		testerMetropolis *= rationalApproximationsMetropolis[i].evaluate(2.).real();
		testerHeatBath *= rationalApproximationsHeatBath[i].evaluate(2.).real();
	}
	double numberFermions = -2*log(testerMetropolis)/log(2.);//We are using the square of the dirac operator, hence we have a factor 2
	if (isOutputProcess()) std::cout << "NFlavorUpdater::The theory has " <<  numberFermions << " nf." << std::endl;
#ifdef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 0.5) < epsilon) std::cout << "NFlavorUpdater::The theory seems SUSY" << std::endl;
#endif
#ifndef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 1) < epsilon) std::cout << "NFlavorUpdater::The theory seems 1 Flavor" << std::endl;
	if (isOutputProcess() && fabs(numberFermions - 2) < epsilon) std::cout << "NFlavorUpdater::The theory seems 2 Flavor" << std::endl;
#endif
	if (isOutputProcess() && fabs(testerMetropolis/testerForce - 1.) > epsilon) {
		std::cout << "NFlavorUpdater::Warning, large mismatch between force and metropolis approximations: " << fabs(testerMetropolis/testerForce - 1.) << std::endl;
	}
	if (isOutputProcess() && fabs(testerMetropolis*testerHeatBath*testerHeatBath - 1.) > epsilon) {
		std::cout << "NFlavorUpdater::Warning, large mismatch between heatbath and metropolis approximations: " << fabs(testerMetropolis*testerHeatBath*testerHeatBath - 1.) << " " << testerMetropolis << " " << testerHeatBath << std::endl;
	}
}

void MultiStepNFlavorUpdater::execute(environment_t& environment) {
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

	//Take the Dirac Operator for Metropolis/Heatbath
	if (diracOperatorMetropolis == 0)  diracOperatorMetropolis = DiracOperator::getInstance(environment.configurations.get<std::string>("MultiStepNFlavorUpdater::dirac_operator_metropolis::dirac_operator"), 1, environment.configurations,  "MultiStepNFlavorUpdater::dirac_operator_metropolis::");
	diracOperatorMetropolis->setLattice(environment.getFermionLattice());

	//Take the Square Dirac Operator for Metropolis/Heatbath
	if (squareDiracOperatorMetropolis == 0) squareDiracOperatorMetropolis = DiracOperator::getInstance(environment.configurations.get<std::string>("MultiStepNFlavorUpdater::dirac_operator_metropolis::dirac_operator"), 2, environment.configurations,  "MultiStepNFlavorUpdater::dirac_operator_metropolis::");
	squareDiracOperatorMetropolis->setLattice(environment.getFermionLattice());

	//Take the Dirac Operator for the force
	if (diracOperatorForce == 0)  diracOperatorForce = DiracOperator::getInstance(environment.configurations.get<std::string>("MultiStepNFlavorUpdater::dirac_operator_metropolis::dirac_operator"), 1, environment.configurations,  "MultiStepNFlavorUpdater::dirac_operator_force::");
	diracOperatorForce->setLattice(environment.getFermionLattice());

	//Take the Square Dirac Operator for the force
	if (squareDiracOperatorForce == 0) squareDiracOperatorForce = DiracOperator::getInstance(environment.configurations.get<std::string>("MultiStepNFlavorUpdater::dirac_operator_metropolis::dirac_operator"), 2, environment.configurations,  "MultiStepNFlavorUpdater::dirac_operator_force::");
	squareDiracOperatorForce->setLattice(environment.getFermionLattice());

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
		j->evaluate(squareDiracOperatorMetropolis, *i, tmp_pseudofermion);
		++j;


		if (i == pseudofermions.begin()) {
			std::string flag = environment.configurations.get<std::string>("check_rational_approximations");

			if (environment.sweep == 0 && environment.iteration == 0 && flag == "true") {
				//Now we test the correctness of rational/heatbath
				long_real_t test = 0.;

				//Now we use the better approximation for the metropolis step
				std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
				//Now we evaluate it with the rational approximation of the inverse
				rational->evaluate(squareDiracOperatorMetropolis,tmp_pseudofermion,*i);
				test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
				if (isOutputProcess()) std::cout << "NFlavoUpdater::Consistency check for the metropolis: " << test - oldPseudoFermionEnergy << std::endl;

				//Now we use the approximation for the force step
				rational = rationalApproximationsForce[0].begin();
				test = 0.;
				//Now we evaluate it with the rational approximation of the inverse
				rational->evaluate(squareDiracOperatorForce,tmp_pseudofermion,*i);
				test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));

				if (isOutputProcess()) std::cout << "NFlavoUpdater::Consistency check for the first level of the force: " << test - oldPseudoFermionEnergy << std::endl;
			}
		}

	}

	//Get the initial energy of momenta
	long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
	//Get the initial energy of the lattice
	long_real_t oldLatticeEnergy = gaugeAction->energy(environment);

	//Set the dirac operator to the new environment
	squareDiracOperatorMetropolis->setLattice(environmentNew.getFermionLattice());
	diracOperatorMetropolis->setLattice(environmentNew.getFermionLattice());
	squareDiracOperatorForce->setLattice(environmentNew.getFermionLattice());
	diracOperatorForce->setLattice(environmentNew.getFermionLattice());

	//Take the fermion action
	if (fermionAction == 0) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		fermionAction = new NFlavorFermionAction*[numberLevels];
		for (int j = 0; j < numberLevels; ++j) {
			fermionAction[j] = new NFlavorFermionAction(squareDiracOperatorForce, diracOperatorForce, rationalApproximationsForce[j]);
			//TODO dirac operator correctly set?
			for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
				//Add the pseudofermions
				fermionAction[j]->addPseudoFermion(&(*i));
			}
			//Set the precision of the inverter
			fermionAction[j]->setForcePrecision(rationalApproximationsForce[j].begin()->getPrecision());
			fermionAction[j]->setForceMaxIterations(environment.configurations.get<unsigned int>("force_inverter_max_steps"));
		}
	}


	//Take the global action
	if (nFlavorAction == 0) nFlavorAction = new NFlavorAction(gaugeAction, fermionAction[0]);//Here we skip the other forces

	//The t-length of a single integration step
	real_t t_length = environment.configurations.get<double>("hmc_t_length");
	std::vector<Force*> forces;
	//The numbers of integration steps
	std::vector<unsigned int> numbers_steps = environment.configurations.get< std::vector<unsigned int> >("number_hmc_steps");
	if (numbers_steps.size() == 1) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 1) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only one time integration only the first level of the force is used!" << std::endl;
		forces.push_back(nFlavorAction);
	} else if (numbers_steps.size() == 2) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 1) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only two time integration only the first level of the force is used!" << std::endl;
		forces.push_back(fermionAction[0]);
		forces.push_back(gaugeAction);
	} else if (numbers_steps.size() == 3) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 2) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only three time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[1]);
		forces.push_back(fermionAction[0]);
		forces.push_back(gaugeAction);
	} else if (numbers_steps.size() == 4) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 3) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only four time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[2]);
		forces.push_back(fermionAction[1]);
		forces.push_back(fermionAction[0]);
		forces.push_back(gaugeAction);
	}
	else {
		if (isOutputProcess()) std::cout << "MultiStepNFlavorHMCUpdater::Warning, NFlavor does not support more than four time integrations!" << std::endl;
		numbers_steps.resize(1);
		forces.push_back(nFlavorAction);
	}

#ifdef DEBUGFORCE
	//Test of the force
	TestForce testForce;
	extended_gauge_lattice_t tmp;
	//Internal calculations needed
	nFlavorAction->updateForce(tmp, environment);
	testForce.genericTestForce(environment, nFlavorAction, tmp[5][2], 5, 2);
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
	diracOperatorMetropolis->setLattice(environmentNew.getFermionLattice());
	squareDiracOperatorMetropolis->setLattice(environmentNew.getFermionLattice());
	//Now we use the better approximation for the metropolis step
	std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		//Now we evaluate it with the rational approximation of the inverse
		rational->evaluate(squareDiracOperatorMetropolis,tmp_pseudofermion,*i);
		newPseudoFermionEnergy += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
		++rational;
	}

	//Global Metropolis Step
	bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy, newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy);

	if (metropolis) {
		environment = environmentNew;
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");
			
			output->push("hmc_history");
			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy));
			output->pop("hmc_history");
			
			output->push("hmc_history");
			output->write("hmc_history", 1);

			output->pop("hmc_history");
		}
	}
	else {
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");
			
			output->push("hmc_history");
			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy));
			output->pop("hmc_history");

			output->push("hmc_history");
			output->write("hmc_history", 0);
			output->pop("hmc_history");

			output->pop("hmc_history");
		}
	}

	delete integrate;
}

void MultiStepNFlavorUpdater::registerParameters(std::map<std::string, Option>& desc) {
	desc["MultiStepNFlavorUpdater::twist"] = Option("MultiStepNFlavorUpdater::twist", 0.0, "set the value of the twist applied to fermions");
	desc["MultiStepNFlavorUpdater::inverter_precision"] = Option("MultiStepNFlavorUpdater::inverter_precision", 1e-11, "set the precision used by the inverter");
	desc["MultiStepNFlavorUpdater::inverter_max_steps"] = Option("MultiStepNFlavorUpdater::inverter_max_steps", 5000, "set the maximum steps used by the inverter");
		
	desc["MultiStepNFlavorUpdater::multigrid"] = Option("MultiStepNFlavorUpdater::multigrid", "false", "Should we use the multigrid inverter? true/false");
	desc["MultiStepNFlavorUpdater::multigrid_basis_dimension"] = Option("MultiStepNFlavorUpdater::multigrid_basis_dimension", 20, "The dimension of the basis for multigrid");
	desc["MultiStepNFlavorUpdater::multigrid_block_size"] = Option("MultiStepNFlavorUpdater::multigrid_block_size", "{4,4,4,4}", "Block size for Multigrid (syntax: {bx,by,bz,bt})");

	desc["MultiStepNFlavorUpdater::sap_block_size"] = Option("MultiStepNFlavorUpdater::sap_block_size", "{4,4,4,4}", "Block size for SAP (syntax: {bx,by,bz,bt})");
	desc["MultiStepNFlavorUpdater::sap_iterations"] = Option("MultiStepNFlavorUpdater::sap_iterations", 5, "The number of sap iterations");
	desc["MultiStepNFlavorUpdater::sap_inverter_precision"] = Option("MultiStepNFlavorUpdater::sap_inverter_precision", 0.00000000001, "The precision of the inner SAP inverter");
	desc["MultiStepNFlavorUpdater::sap_inverter_max_steps"] = Option("MultiStepNFlavorUpdater::sap_inverter_max_steps", 50, "The maximum number of steps for the inner SAP inverter");
	desc["MultiStepNFlavorUpdater::gmres_inverter_precision"] = Option("MultiStepNFlavorUpdater::gmres_inverter_precision", 0.00000000001, "The precision of the GMRES inverter used to initialize the multigrid basis");
	desc["MultiStepNFlavorUpdater::gmres_inverter_max_steps"] = Option("MultiStepNFlavorUpdater::gmres_inverter_max_steps", 100, "The maximum number of steps for the GMRES inverter used to initialize the multigrid basis");
	DiracOperator::registerParameters(desc, "MultiStepNFlavorUpdater::dirac_operator_metropolis::");
	DiracOperator::registerParameters(desc, "MultiStepNFlavorUpdater::dirac_operator_force::");
}

} /* namespace Update */
