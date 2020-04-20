#include "SingletOperators.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/StoutSmearing.h"
#include "utils/Gamma.h"
#include "dirac_operators/BlockDiracWilsonOperator.h"
#include "dirac_operators/ComplementBlockDiracOperator.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "multigrid/MultiGridSolver.h"
#include "inverters/GMRESR.h"
#include "dirac_operators/GammaOperators.h"
#include "dirac_operators/HoppingOperator.h"

namespace Update {

SingletOperators::SingletOperators() : LatticeSweep(), MultiGridStochasticEstimator(), diracOperator(0), squareDiracOperator(0), inverter(0), gamma() { }

SingletOperators::SingletOperators(const SingletOperators& toCopy) : LatticeSweep(toCopy), MultiGridStochasticEstimator(toCopy), diracOperator(0), squareDiracOperator(0), inverter(0), gamma() { }

SingletOperators::~SingletOperators() {
	if (diracOperator) delete diracOperator;
}

void SingletOperators::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;//TODO: only vector operations?

	//unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");
	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}
	if (squareDiracOperator == 0) {
		squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	}

	diracOperator->setLattice(environment.getFermionLattice());
	diracOperator->setGamma5(false);
	squareDiracOperator->setLattice(environment.getFermionLattice());

	//Here we construct the SAP preconditioner
	BlockDiracOperator* blackBlockDiracOperator = 0;
	BlockDiracOperator* redBlockDiracOperator = 0;
	
	if (environment.configurations.get<std::string>("SingletOperators::multigrid") == "true") {
		blackBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Black);
		blackBlockDiracOperator->setLattice(environment.getFermionLattice());
		blackBlockDiracOperator->setGamma5(false);
		blackBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("SingletOperators::sap_block_size"));
		
		redBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Red);
		redBlockDiracOperator->setLattice(environment.getFermionLattice());
		redBlockDiracOperator->setGamma5(false);
		redBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("SingletOperators::sap_block_size"));

		MultiGridSolver* multiGridSolver = new MultiGridSolver(environment.configurations.get< unsigned int >("SingletOperators::multigrid_basis_dimension"), environment.configurations.get< std::vector<unsigned int> >("SingletOperators::multigrid_block_size"), blackBlockDiracOperator, redBlockDiracOperator);
		multiGridSolver->setSAPIterations(environment.configurations.get<unsigned int>("SingletOperators::sap_iterations"));
		multiGridSolver->setSAPMaxSteps(environment.configurations.get<unsigned int>("SingletOperators::sap_inverter_max_steps"));
		multiGridSolver->setSAPPrecision(environment.configurations.get<real_t>("SingletOperators::sap_inverter_precision"));
		multiGridSolver->setGMRESIterations(environment.configurations.get<unsigned int>("SingletOperators::gmres_inverter_max_steps"));
		multiGridSolver->setGMRESPrecision(environment.configurations.get<real_t>("SingletOperators::gmres_inverter_precision"));

		multiGridSolver->initializeBasis(diracOperator);

		inverter = multiGridSolver;

		if (isOutputProcess()) std::cout << "SingletOperators::Using multigrid inverter and SAP preconditioning ..." << std::endl;
	}
	else {
		GMRESR* gmresr = new GMRESR();
		inverter = gmresr;
		

		if (isOutputProcess()) std::cout << "SingletOperators::Without using multigrid ..." << std::endl;
	}

	inverter->setPrecision(environment.configurations.get<double>("SingletOperators::inverter_precision"));
	inverter->setMaximumSteps(environment.configurations.get<unsigned int>("SingletOperators::inverter_max_steps"));

	int inversionSteps = 0;

	
	unsigned int numberRandomSources = environment.configurations.get<unsigned int>("SingletOperators::number_stochastic_estimators");
	extended_dirac_vector_t* randomSources = new extended_dirac_vector_t[numberRandomSources];
	extended_dirac_vector_t* inverseRandomSources = new extended_dirac_vector_t[numberRandomSources];

	MultiGridSolver* mgs = dynamic_cast<MultiGridSolver*>(inverter);
	if (mgs && environment.configurations.get<std::string>("SingletOperators::use_multigrid_diluition") == "true") {
		this->setBlockBasis(mgs->getBasis());
		for (unsigned int i = 0; i < numberRandomSources; ++i) {
			if (i % 3 != 0) {
				this->getMultigridVectors(diracOperator, randomSources[i], inverseRandomSources[i]);
			}
			else {
				this->getResidualVectors(inverter, diracOperator, randomSources[i], inverseRandomSources[i]);
			}
		}
	} else {
		extended_dirac_vector_t tmp;
		if (isOutputProcess()) std::cout << "SingletOperators::Disconnected contributions measured without multigrid" << std::endl;
		for (unsigned int i = 0; i < numberRandomSources; ++i) {
			this->generateRandomNoise(randomSources[i]);
			inverter->solve(diracOperator, randomSources[i], inverseRandomSources[i]);
		}
	}

	HoppingOperator H(diracOperator);
	GammaOperators gammaOperators;
	extended_dirac_vector_t tmp;

	for (int hopping = 0; hopping < 9; ++hopping) {
		std::complex<long_real_t> disconnectedResults[numberRandomSources][16];

		for (unsigned int source = 0; source < numberRandomSources; ++source) {
			for (int  v = 0; v < 16; ++v) {
				gammaOperators.multiply(tmp,inverseRandomSources[source], v);
				disconnectedResults[source][v] = AlgebraUtils::dot(randomSources[source],tmp);
					
				//Correct normalization
				long_real_t factor = 2*diracOperator->getKappa();
				disconnectedResults[source][v] = disconnectedResults[source][v]*factor;
			}
		}
	
		if (environment.measurement && isOutputProcess() && hopping == 0) {
			GlobalOutput* output = GlobalOutput::getInstance();
		
			output->push("disconnected_operators");
			for (int v = 0; v < 16; ++v) {
				output->push("disconnected_operators");
				std::cout << "Disconnected operator " << v << ": ";
	
				for (unsigned int i = 0; i < numberRandomSources; ++i) {
					if (i != numberRandomSources - 1) std::cout << disconnectedResults[i][v] << ", ";
					else std::cout << disconnectedResults[i][v] << std::endl;
				
					output->push("disconnected_operators");
					output->write("disconnected_operators", real(disconnectedResults[i][v]));
					output->write("disconnected_operators", imag(disconnectedResults[i][v]));
					output->pop("disconnected_operators");
				}
				output->pop("disconnected_operators");
			}
			output->pop("disconnected_operators");
		}

		if (environment.measurement && isOutputProcess() && hopping != 0) {
			GlobalOutput* output = GlobalOutput::getInstance();
			std::string name = "disconnected_operators_hopping_";
			name += toString(hopping)+"_dirac_inverse";

			output->push(name);
			for (int v = 0; v < 16; ++v) {
				output->push(name);
				std::cout << "Disconnected operator " << v << " hopping " << hopping << ": ";

				for (unsigned int i = 0; i < numberRandomSources; ++i) {
					if (i != numberRandomSources - 1) std::cout << disconnectedResults[i][v] << ", ";
					else std::cout << disconnectedResults[i][v] << std::endl;
				
					output->push(name);
					output->write(name, real(disconnectedResults[i][v]));
					output->write(name, imag(disconnectedResults[i][v]));
					output->pop(name);
				}

				output->pop(name);
			}
			output->pop(name);
		}
			
		for (unsigned int i = 0; i < numberRandomSources; ++i) {
			H.apply(tmp, inverseRandomSources[i]);
			inverseRandomSources[i] = tmp;
		}
	}


	//Here we evaluate the trace of the hopping parameter expansion terms
	for (int hopping = 1; hopping < 8; ++hopping) {
			
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			std::string name = "disconnected_operators_hopping_";
			name += toString(hopping);

			output->push(name);
		}
	}

	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);

	extended_dirac_vector_t source, hopping_source, previous;
	unsigned int hopping_stochastic_estimators = environment.configurations.get<unsigned int>("SingletOperators::number_stochastic_estimators_hopping_terms");
				

	for (unsigned int i = 0; i < hopping_stochastic_estimators; ++i) {
		this->generateRandomNoise(source);
		previous = source;
		for (int hopping = 1; hopping < 8; ++hopping) {
			H.apply(hopping_source, previous);
			previous = hopping_source;

			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				std::string name = "disconnected_operators_hopping_";
				name += toString(hopping);

				output->push(name);
			}

			for (int  v = 0; v < 16; ++v) {
				gammaOperators.multiply(tmp,hopping_source, v);
				std::complex<long_real_t> result = AlgebraUtils::dot(source,tmp);

				//Correct normalization
				long_real_t factor = 2*diracOperator->getKappa();
				result = factor*result;

				if (environment.measurement && isOutputProcess()) {
					GlobalOutput* output = GlobalOutput::getInstance();
					std::string name = "disconnected_operators_hopping_";
					name += toString(hopping);

					output->push(name);

					output->write(name, real(result));
					output->write(name, imag(result));

					output->pop(name);
				}
			}

			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();
				std::string name = "disconnected_operators_hopping_";
				name += toString(hopping);

				output->pop(name);
			}
			
		}
	}

	for (int hopping = 1; hopping < 8; ++hopping) {
		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			std::string name = "disconnected_operators_hopping_";
			name += toString(hopping);
			output->pop(name);
		}
	}
		
	delete[] randomSources;
	delete[] inverseRandomSources;
	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	if (isOutputProcess()) std::cout << "SingletOperators::Measure of hopping parameter expansion done in: " << (elapsed) << " s."<< std::endl;
	
	if (redBlockDiracOperator) delete redBlockDiracOperator;
	if (blackBlockDiracOperator) delete blackBlockDiracOperator;
	if (inverter) delete inverter;
}

void SingletOperators::registerParameters(po::options_description& desc) {
	desc.add_options()
		("SingletOperators::inverter_precision", po::value<double>()->default_value(0.0000000001), "set the precision used by the inverter")
		("SingletOperators::inverter_max_steps", po::value<unsigned int>()->default_value(5000), "set the maximum steps used by the inverter")
		
		("SingletOperators::multigrid", po::value<std::string>()->default_value("true"), "Should we use the multigrid inverter? true/false")
		("SingletOperators::multigrid_basis_dimension", po::value<unsigned int>()->default_value(20), "The dimension of the basis for multigrid")
		("SingletOperators::multigrid_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for Multigrid (syntax: {bx,by,bz,bt})")

		("SingletOperators::sap_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for SAP (syntax: {bx,by,bz,bt})")
		("SingletOperators::sap_iterations", po::value<unsigned int>()->default_value(5), "The number of sap iterations")
		("SingletOperators::sap_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the inner SAP inverter")
		("SingletOperators::sap_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the inner SAP inverter")
		("SingletOperators::gmres_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the GMRES inverter used to initialize the multigrid basis")
		("SingletOperators::gmres_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the GMRES inverter used to initialize the multigrid basis")
		
		("SingletOperators::number_stochastic_estimators", po::value<unsigned int>()->default_value(13), "Number of stochastic estimators for the disconnected part")
		("SingletOperators::number_stochastic_estimators_hopping_terms", po::value<unsigned int>()->default_value(2500), "Number of stochastic estimators for the disconnected part in the hopping parameter expansion")
		("SingletOperators::use_multigrid_diluition", po::value<std::string>()->default_value("true"), "Should we use the multigrid diluition for the measure of disconnected contribution? true/false")
		;
}

} /* namespace Update */
