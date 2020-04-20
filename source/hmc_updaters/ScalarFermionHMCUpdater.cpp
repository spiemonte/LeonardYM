#include "ScalarFermionHMCUpdater.h"
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

ScalarFermionHMCUpdater::ScalarFermionHMCUpdater() : LatticeSweep(), nFlavorQCDAction(0), gaugeAction(0), fermionAction(0), scalarAction(0), squareDiracOperator(0), diracOperator(0), multishiftSolver(0) { }

ScalarFermionHMCUpdater::ScalarFermionHMCUpdater(const ScalarFermionHMCUpdater& toCopy) : LatticeSweep(toCopy), nFlavorQCDAction(0), gaugeAction(0), fermionAction(0), scalarAction(0), squareDiracOperator(0), diracOperator(0), multishiftSolver(0) { }

ScalarFermionHMCUpdater::~ScalarFermionHMCUpdater() {
	if (scalarAction != 0) delete scalarAction;
	if (nFlavorQCDAction != 0) delete nFlavorQCDAction;
	if (multishiftSolver != 0) delete multishiftSolver;
}

void ScalarFermionHMCUpdater::initializeApproximations(environment_t& environment) {
	double twist = 0;
	try {
		twist = environment.configurations.get<double>("ScalarFermionHMCUpdater::twist");
	} catch(NotFoundOption) {
		
	}
	if (isOutputProcess()) std::cout << "ScalarFermionHMCUpdater::Using twist " << twist << std::endl;	

	if (multishiftSolver == 0) {
		multishiftSolver = MultishiftSolver::getInstance("minimal_residual");
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
	try {
		for (int i = 1; i <= numberLevels; ++i) {
			level_precisions[i-1] = environment.configurations.get<double>(std::string("force_inverter_precision_level_")+toString(i));
		}
	} catch (NotFoundOption& ex) {
		for (int i = 1; i <= numberLevels; ++i) {
			level_precisions[i-1] = environment.configurations.get<double>("force_inverter_precision");
		}
		if (isOutputProcess()) std::cout << "ScalarFermionHMCUpdater::Warning, a single precision is provided for all the level of the force!" << std::endl;
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

void ScalarFermionHMCUpdater::checkTheory(const environment_t& environment) const {
	double testerForce = 1., testerMetropolis = 1., testerHeatBath = 1.;
	double epsilon = 0.0001;
	int numberPseudofermions = environment.configurations.get< unsigned int >("number_pseudofermions");
	for (int i = 0; i < numberPseudofermions; ++i) {
		testerForce *= rationalApproximationsForce[0][i].evaluate(2.).real();
		testerMetropolis *= rationalApproximationsMetropolis[i].evaluate(2.).real();
		testerHeatBath *= rationalApproximationsHeatBath[i].evaluate(2.).real();
	}
	double numberFermions = -2*log(testerMetropolis)/log(2.);//We are using the square of the dirac operator, hence we have a factor 2
	if (isOutputProcess()) std::cout << "ScalarFermionHCMUpdater::The theory has " <<  numberFermions << " fermions." << std::endl;
#ifdef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 0.5) < epsilon) std::cout << "ScalarFermionHMCUpdater::The fermion sector of the theory seems SUSY" << std::endl;
#endif
#ifndef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 1) < epsilon) std::cout << "ScalarFermionHMCUpdater::The fermion sector of the theory seems 1 FlavorQCD" << std::endl;
	if (isOutputProcess() && fabs(numberFermions - 2) < epsilon) std::cout << "ScalarFermionHMCUpdater::The fermion sector of the theory seems 2 FlavorQCD" << std::endl;
#endif
	if (isOutputProcess() && fabs(testerMetropolis/testerForce - 1.) > epsilon) {
		std::cout << "ScalarFermionUpdater::Warning, large mismatch between force and metropolis approximations: " << fabs(testerMetropolis/testerForce - 1.) << std::endl;
	}
	if (isOutputProcess() && fabs(testerMetropolis*testerHeatBath*testerHeatBath - 1.) > epsilon) {
		std::cout << "ScalarFermionUpdater::Warning, large mismatch between heatbath and metropolis approximations: " << fabs(testerMetropolis*testerHeatBath*testerHeatBath - 1.) << " " << testerMetropolis << " " << testerHeatBath << std::endl;
	}
}

void ScalarFermionHMCUpdater::initializeScalarAction(environment_t& environment) {
	unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");
	unsigned int nf = environment.configurations.get<unsigned int>("fundamental_nf_scalars");

	if (scalarAction == 0 ) {
		scalarAction = new MultiScalarAction();
		if (aNf != 0) {
                	double mu = environment.configurations.get<double>("ScalarFermionHMCUpdater::adjoint_scalar_mass");
                	double lambda = environment.configurations.get<double>("ScalarFermionHMCUpdater::adjoint_quartic_coupling");
			double lambda_8 = environment.configurations.get<double>("ScalarFermionHMCUpdater::adjoint_adjoint_quartic_coupling");

                	AdjointScalarAction* adjointScalarAction = new AdjointScalarAction(mu, lambda, lambda_8);

			scalarAction->addAction(adjointScalarAction);
		}
		if (nf != 0) {
			double mu = environment.configurations.get<double>("ScalarFermionHMCUpdater::fundamental_scalar_mass");
                        double lambda = environment.configurations.get<double>("ScalarFermionHMCUpdater::fundamental_quartic_coupling");
			double lambda_8 = environment.configurations.get<double>("ScalarFermionHMCUpdater::fundamental_adjoint_quartic_coupling");

                        FundamentalScalarAction* fundamentalScalarAction = new FundamentalScalarAction(mu, lambda, lambda_8);

                        scalarAction->addAction(fundamentalScalarAction);
		}
        }

        if (environment.adjoint_scalar_fields.size() != aNf) {
                if (isOutputProcess()) std::cout << "ScalarFermionHMCUpdater::Error, adjoint scalar fields not initialized, set " << aNf << ", but now " << environment.adjoint_scalar_fields.size() << std::endl;
		exit(97);
        }
	if (environment.fundamental_scalar_fields.size() != nf) {
                if (isOutputProcess()) std::cout << "ScalarFermionHMCUpdater::Error, fundamental scalar fields not initialized, set " << nf << ", but now " << environment.fundamental_scalar_fields.size() << std::endl;
                exit(97);
        }
}

void ScalarFermionHMCUpdater::execute(environment_t& environment) {
	//First we initialize the approximations
	this->initializeApproximations(environment);

	//Then we initialize the scalar action
	this->initializeScalarAction(environment);

	//We check the theory that is simulated
	if (environment.iteration == 0 && environment.sweep == 0) {
		this->checkTheory(environment);
	}

	//Initialize the momenta
	this->randomMomenta(momenta);
	//Copy the environment
	environmentNew = environment;

	//Take the Dirac Operator
	if (diracOperator == 0)  diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());

	//Take the Square Dirac Operator
	if (squareDiracOperator == 0) squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
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
		j->evaluate(squareDiracOperator, *i, tmp_pseudofermion);
		++j;


		if (i == pseudofermions.begin()) {
			try {
				std::string flag = environment.configurations.get<std::string>("check_rational_approximations");

				if (environment.sweep == 0 && environment.iteration == 0 && flag == "true") {
					//Now we test the correctness of rational/heatbath
					long_real_t test = 0.;

					//Now we use the better approximation for the metropolis step
					std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
					//Now we evaluate it with the rational approximation of the inverse
					rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
					test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
					if (isOutputProcess()) std::cout << "NFlavoQCDUpdater::Consistency check for the metropolis: " << test - oldPseudoFermionEnergy << std::endl;

					//Now we use the approximation for the force step
					rational = rationalApproximationsForce[0].begin();
					test = 0.;
					//Now we evaluate it with the rational approximation of the inverse
					rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
					test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));

					if (isOutputProcess()) std::cout << "NFlavoQCDUpdater::Consistency check for the first level of the force: " << test - oldPseudoFermionEnergy << std::endl;
				}

			} catch (NotFoundOption& ex) {
				if (isOutputProcess() && environment.sweep == 0 && environment.iteration == 0 && environment.measurement) std::cout << "NFlavorQCDUpdater::No consistency check of metropolis/force approximations!" << std::endl;
			}
		}

	}

	//Get the initial energy of momenta
	long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
	//Get the initial energy of the lattice
	long_real_t oldLatticeEnergy = gaugeAction->energy(environment);
	//Get the action of the scalars
	long_real_t oldScalarEnergy = scalarAction->energy(environment);

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
			fermionAction[j]->setForcePrecision(rationalApproximationsForce[j].begin()->getPrecision());
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

	if (numbers_steps.size() == 2) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 1) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only two time integration only the first level of the force is used!" << std::endl;
		forces.push_back(nFlavorQCDAction);
		forces.push_back(scalarAction);
	} else if (numbers_steps.size() == 3) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 2) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only three time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[0]);
		forces.push_back(scalarAction);
		forces.push_back(gaugeAction);
	} else if (numbers_steps.size() == 4) {
		int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
		if (isOutputProcess() && numberLevels != 3) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only four time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[1]);
		forces.push_back(fermionAction[0]);
		forces.push_back(scalarAction);
		forces.push_back(gaugeAction);
	} else if (numbers_steps.size() == 5) {
                int numberLevels = environment.configurations.get< unsigned int >("number_force_levels");
                if (isOutputProcess() && numberLevels != 3) std::cout << "MultiStepNFlavorHMCUpdater::Warning, with only four time integration only the first two levels of the force is used!" << std::endl;
		forces.push_back(fermionAction[2]);
                forces.push_back(fermionAction[1]);
                forces.push_back(fermionAction[0]);
                forces.push_back(scalarAction);
                forces.push_back(gaugeAction);
        }
	else {
		if (isOutputProcess()) std::cout << "FermionScalarUpdater::Warning, NFlavor does not support more than four time or only one integration step!" << std::endl;
		numbers_steps.resize(1);
		forces.push_back(nFlavorQCDAction);
		forces.push_back(scalarAction);
	}

#ifdef DEBUGFORCE
	//Test of the force
	TestForce testForce;
	extended_gauge_lattice_t tmp;
	//Internal calculations needed
	scalarAction->updateForce(tmp, environment);
	testForce.genericTestForce(environment, scalarAction, tmp[5][2], 5, 2);
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
	//Get the final scalar energy of the lattice
	long_real_t newScalarEnergy = scalarAction->energy(environmentNew);
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
	bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy + oldScalarEnergy, newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy + newScalarEnergy);

	if (metropolis) {
		environment = environmentNew;
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");

			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy + oldScalarEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy + newScalarEnergy));
			output->write("hmc_history", 1);

			output->pop("hmc_history");
		}
	}
	else {
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("hmc_history");

			output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy + oldScalarEnergy) + (newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy + newScalarEnergy));
			output->write("hmc_history", 0);

			output->pop("hmc_history");
		}
	}

	delete integrate;
}

void ScalarFermionHMCUpdater::registerParameters(po::options_description& desc) {
	static bool single = true;
	if (single) desc.add_options()
		("ScalarFermionHMCUpdater::adjoint_scalar_mass", po::value<double>()->default_value(0.0), "set the value of the adjoint scalar mass m*m")
		("ScalarFermionHMCUpdater::adjoint_quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the adjoint quartic coupling lambda")
		("ScalarFermionHMCUpdater::adjoint_adjoint_quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the adjoint quartic coupling lambda")
		("ScalarFermionHMCUpdater::fundamental_scalar_mass", po::value<double>()->default_value(0.0), "set the value of the adjoint scalar mass m*m")
                ("ScalarFermionHMCUpdater::fundamental_quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the adjoint quartic coupling lambda")
		("ScalarFermionHMCUpdater::fundamental_adjoint_quartic_coupling", po::value<double>()->default_value(0.0), "set the value of the adjoint quartic coupling lambda")
		("ScalarFermionHMCUpdater::twist", po::value<double>()->default_value(0.0), "set the value of the twist applied to fermions")
		("ScalarFermionHMCUpdater::inverter_precision", po::value<double>()->default_value(0.0000000001), "set the precision used by the inverter")
		("ScalarFermionHMCUpdater::inverter_max_steps", po::value<unsigned int>()->default_value(5000), "set the maximum steps used by the inverter")
		;
	single = false;
}

} /* namespace Update */
