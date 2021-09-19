#include "HiggsGaugeHMCUpdater.h"
#include "hmc_integrators/Integrate.h"
#include "actions/GaugeAction.h"
#include "hmc_forces/TestForce.h"
#include "io/GlobalOutput.h"
#include <string>

//#define REVERSIBILITY_CHECK
//#define DEBUGFORCE

namespace Update {

	HiggsGaugeHMCUpdater::HiggsGaugeHMCUpdater() : LatticeSweep(), HMCUpdater(), scalarAction(0) { }

	HiggsGaugeHMCUpdater::HiggsGaugeHMCUpdater(const HiggsGaugeHMCUpdater& toCopy) : LatticeSweep(toCopy), HMCUpdater(toCopy), scalarAction(0) { }

	HiggsGaugeHMCUpdater::~HiggsGaugeHMCUpdater() {
		if (scalarAction != 0) delete scalarAction;
	}

	void HiggsGaugeHMCUpdater::initializeScalarAction(environment_t& environment) {
		unsigned int aNf = environment.configurations.get<unsigned int>("adjoint_nf_scalars");
		unsigned int nf = environment.configurations.get<unsigned int>("fundamental_nf_scalars");

		if (scalarAction == 0) {
			double adjoint_scalar_mass = environment.configurations.get<double>("HiggsGaugeHMCUpdater::adjoint_scalar_mass");
			double fundamental_scalar_mass = environment.configurations.get<double>("HiggsGaugeHMCUpdater::fundamental_scalar_mass");
			double lambda_adjoint = environment.configurations.get<double>("HiggsGaugeHMCUpdater::adjoint_quartic_coupling");
			double lambda_fundamental = environment.configurations.get<double>("HiggsGaugeHMCUpdater::fundamental_quartic_coupling");
			double lambda_8 = environment.configurations.get<double>("HiggsGaugeHMCUpdater::lambda_8");
			double lambda_mixed = environment.configurations.get<double>("HiggsGaugeHMCUpdater::lambda_mixed");

			scalarAction = new ScalarAction(adjoint_scalar_mass, fundamental_scalar_mass, lambda_fundamental, lambda_adjoint, lambda_mixed, lambda_8);
		}

		if (environment.adjoint_scalar_fields.size() != aNf) {
			if (isOutputProcess()) std::cout << "HiggsGaugeHMCUpdater::Error, adjoint scalar fields not initialized, set " << aNf << ", but now " << environment.adjoint_scalar_fields.size() << std::endl;
			exit(97);
		}
		if (environment.fundamental_scalar_fields.size() != nf) {
			if (isOutputProcess()) std::cout << "HiggsGaugeHMCUpdater::Error, fundamental scalar fields not initialized, set " << nf << ", but now " << environment.fundamental_scalar_fields.size() << std::endl;
			exit(97);
		}
	}

	void HiggsGaugeHMCUpdater::execute(environment_t& environment) {
		//Then we initialize the scalar action
		this->initializeScalarAction(environment);

		//Initialize the momenta
		this->randomMomenta(momenta);

		//Copy the lattice and the environment
		environmentNew = environment;

		real_t beta = environment.configurations.get<double>("beta");

		//Get the gauge action
		GaugeAction* gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),beta);

		//Get the initial energy of momenta
		long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
		//Get the initial energy of the lattice
		long_real_t oldLatticeEnergy = gaugeAction->energy(environment);
		//Get the action of the scalars
		long_real_t oldScalarEnergy = scalarAction->energy(environment);

		//The t-length of a single integration step
		real_t t_length = environment.configurations.get<double>("hmc_t_length");
		//The numbers of integration steps
		std::vector<unsigned int> numbers_steps = environment.configurations.get< std::vector<unsigned int> >("number_hmc_steps");
		if (numbers_steps.size() != 2) {
			if (isOutputProcess()) std::cout << "HiggsGaugeHMCUpdater::Warning, higgs-gauge does support only two time integrations!" << std::endl;
			numbers_steps.resize(2);
			if (numbers_steps[0] == 0) numbers_steps[0] = 1;
			if (numbers_steps[1] == 0) numbers_steps[1] = 1;
		}

		Integrate* integrate = Integrate::getInstance(environment.configurations.get<std::string>("name_integrator"));

		//The vector of the force for the integrator
		std::vector<Force*> force;
		force.push_back(scalarAction);
		force.push_back(gaugeAction);


#ifdef DEBUGFORCE
        //Test of the force
		TestForce testForce;
		extended_gauge_lattice_t tmp;
        //Internal calculations needed
		scalarAction->updateForce(tmp, environment);
		testForce.genericTestForce(environment, scalarAction, tmp[5][2], 5, 2);
#endif

		//Integrate numerically the equation of motion
		integrate->integrate(environmentNew, momenta, force, numbers_steps, t_length);
#ifdef REVERSIBILITY_CHECK
		integrate->integrate(environmentNew, momenta, force, numbers_steps, -t_length);
#endif

		//Get the final energy of momenta
		long_real_t newMomentaEnergy = this->momentaEnergy(momenta);
		//Get the final energy of the lattice
		long_real_t newLatticeEnergy = gaugeAction->energy(environmentNew);
		//Get the final scalar energy of the lattice
		long_real_t newScalarEnergy = scalarAction->energy(environmentNew);

#ifdef REVERSIBILITY_CHECK
		if (isOutputProcess()) {
			std::cout << "HiggsGaugeHMCUpdater::momenta energy difference: " << newMomentaEnergy - oldMomentaEnergy << std::endl;
			std::cout << "HiggsGaugeHMCUpdater::gauge energy difference: " << newLatticeEnergy - oldLatticeEnergy << std::endl;
			std::cout << "HiggsGaugeHMCUpdater::scalar energy difference: " << newScalarEnergy - oldScalarEnergy << std::endl;
		}
#endif

		//Global Metropolis Step
		bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldScalarEnergy, newMomentaEnergy + newLatticeEnergy + newScalarEnergy);

		if (metropolis) {
			environment = environmentNew;
			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();

				output->push("hmc_history");

				output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldScalarEnergy) + (newMomentaEnergy + newLatticeEnergy + newScalarEnergy));
				output->write("hmc_history", 1);

				output->pop("hmc_history");
			}
		}
		else {
			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();

				output->push("hmc_history");

				output->write("hmc_history", - (oldMomentaEnergy + oldLatticeEnergy + oldScalarEnergy) + (newMomentaEnergy + newLatticeEnergy + newScalarEnergy));
				output->write("hmc_history", 0);

				output->pop("hmc_history");
			}
		}

		delete integrate;
		delete gaugeAction;
	}

	void HiggsGaugeHMCUpdater::registerParameters(std::map<std::string, Option>& desc) {
		desc["HiggsGaugeHMCUpdater::adjoint_scalar_mass"] = Option("HiggsGaugeHMCUpdater::adjoint_scalar_mass", 0.0, "set the value of the adjoint scalar mass m*m");
        desc["HiggsGaugeHMCUpdater::fundamental_scalar_mass"] = Option("HiggsGaugeHMCUpdater::fundamental_scalar_mass", 0.0, "set the value of the fundamental scalar mass m*m");

        desc["HiggsGaugeHMCUpdater::adjoint_quartic_coupling"] = Option("HiggsGaugeHMCUpdater::adjoint_quartic_coupling", 0.0, "set the value of the adjoint quartic coupling lambda");
        desc["HiggsGaugeHMCUpdater::fundamental_quartic_coupling"] = Option("HiggsGaugeHMCUpdater::fundamental_quartic_coupling", 0.0, "set the value of the fundamental quartic coupling lambda");

        desc["HiggsGaugeHMCUpdater::lambda_8"] = Option("HiggsGaugeHMCUpdater::lambda_8", 0.0, "set the value of the adjoint quartic coupling lambda");
        desc["HiggsGaugeHMCUpdater::lambda_mixed"] = Option("HiggsGaugeHMCUpdater::lambda_mixed", 0.0, "set the value of the mixed quartic coupling lambda");
	}

} /* namespace Update */
