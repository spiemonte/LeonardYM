#include "NPRVertex.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/StoutSmearing.h"
#include "utils/Gamma.h"
#include "dirac_operators/BlockDiracOperator.h"
#include "dirac_operators/ComplementBlockDiracOperator.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "multigrid/MultiGridSolver.h"
#include "inverters/GMRESR.h"
#include "inverters/PreconditionedBiCGStab.h"
#include "dirac_operators/GammaOperators.h"
#include "dirac_operators/HoppingOperator.h"

namespace Update {

NPRVertex::NPRVertex() : LatticeSweep(), MultiGridStochasticEstimator(), diracOperator(0), inverter(0), gamma() { }

NPRVertex::NPRVertex(const NPRVertex& toCopy) : LatticeSweep(toCopy), MultiGridStochasticEstimator(toCopy), diracOperator(0), inverter(0), gamma() { }

NPRVertex::~NPRVertex() {
	if (diracOperator) delete diracOperator;
}

void NPRVertex::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;//TODO: only vector operations?

	//unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");
	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	diracOperator->setLattice(environment.getFermionLattice());
	diracOperator->setGamma5(false);

	//Here we construct the SAP preconditioner
	BlockDiracOperator* blackBlockDiracOperator = 0;
	BlockDiracOperator* redBlockDiracOperator = 0;
	
	if (environment.configurations.get<std::string>("NPRVertex::multigrid") == "true") {
		blackBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Black);
		blackBlockDiracOperator->setLattice(environment.getFermionLattice());
		blackBlockDiracOperator->setGamma5(false);
		blackBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("NPRVertex::sap_block_size"));
		
		redBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Red);
		redBlockDiracOperator->setLattice(environment.getFermionLattice());
		redBlockDiracOperator->setGamma5(false);
		redBlockDiracOperator->setBlockSize(environment.configurations.get< std::vector<unsigned int> >("NPRVertex::sap_block_size"));

		MultiGridSolver* multiGridSolver = new MultiGridSolver(environment.configurations.get< unsigned int >("NPRVertex::multigrid_basis_dimension"), environment.configurations.get< std::vector<unsigned int> >("NPRVertex::multigrid_block_size"), blackBlockDiracOperator, redBlockDiracOperator);
		multiGridSolver->setSAPIterations(environment.configurations.get<unsigned int>("NPRVertex::sap_iterations"));
		multiGridSolver->setSAPMaxSteps(environment.configurations.get<unsigned int>("NPRVertex::sap_inverter_max_steps"));
		multiGridSolver->setSAPPrecision(environment.configurations.get<real_t>("NPRVertex::sap_inverter_precision"));
		multiGridSolver->setGMRESIterations(environment.configurations.get<unsigned int>("NPRVertex::gmres_inverter_max_steps"));
		multiGridSolver->setGMRESPrecision(environment.configurations.get<real_t>("NPRVertex::gmres_inverter_precision"));

		multiGridSolver->initializeBasis(diracOperator);

		inverter = multiGridSolver;

		if (isOutputProcess()) std::cout << "NPRVertex::Using multigrid inverter and SAP preconditioning ..." << std::endl;
	}
	else {
		PreconditionedBiCGStab* pbicg = new PreconditionedBiCGStab();
		inverter = pbicg;
		

		if (isOutputProcess()) std::cout << "NPRVertex::Without using multigrid ..." << std::endl;
	}

	inverter->setPrecision(environment.configurations.get<double>("NPRVertex::inverter_precision"));
	inverter->setMaximumSteps(environment.configurations.get<unsigned int>("NPRVertex::inverter_max_steps"));

	std::complex<long_real_t> vertex[16][4*diracVectorLength][4*diracVectorLength];

	int inversionSteps = 0;

	std::vector<real_t> momentum = environment.configurations.get< std::vector<real_t> >("NPRVertex::momentum");
	if (isOutputProcess()) std::cout << "NPRVertex::Momentum p = {" << momentum[0] << ", " << momentum[1] << ", " << momentum[2] << ", " << momentum[3] << "}" << std::endl;
	momentum[0] = momentum[0]*2.*PI/Layout::glob_x;
	momentum[1] = momentum[1]*2.*PI/Layout::glob_y;
	momentum[2] = momentum[2]*2.*PI/Layout::glob_z;
	momentum[3] = momentum[3]*2.*PI/Layout::glob_t;
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateMomentumSource(source[c*4 + alpha], momentum, alpha, c);
			//this->generateSource(source, alpha, c);
			inverter->solve(diracOperator, source[c*4 + alpha], inverse_source[c*4 + alpha]);
			
			extended_dirac_vector_t test;
			diracOperator->multiply(test,inverse_source[c*4 + alpha]);
			long_real_t dtest = AlgebraUtils::differenceNorm(test, source[c*4 + alpha]);
			if (isOutputProcess()) std::cout << "NPRVertex::Convergence test of the inverter : " << dtest << std::endl;

			inversionSteps += inverter->getLastSteps();
		}
	}
	
	if (isOutputProcess()) std::cout << "NPRVertex::Vertex computed with " << inversionSteps << " inversion steps" << std::endl;
	

	long_real_t factor = 2.*diracOperator->getKappa();

	//First we compute the propagator
	std::complex<long_real_t> propagator[4*diracVectorLength][4*diracVectorLength];
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			propagator[i][j] = factor*AlgebraUtils::dot(source[i], inverse_source[j])/static_cast<long_real_t>(Layout::globalVolume);
		}
	}

	//Now we compute the vertex
	//Ugly initialization to zero
#ifdef MULTITHREADING
	std::complex<long_real_t> resultVertex[16][4*diracVectorLength][4*diracVectorLength][omp_get_max_threads()];
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			for (int k = 0; k < omp_get_max_threads(); ++k) {
				for (int v = 0; v < 16; ++v) {
					resultVertex[v][i][j][k] = 0.;
				}
			}
		}
	}
#endif
#ifndef MULTITHREADING
	std::complex<long_real_t> resultVertex[16][4*diracVectorLength][4*diracVectorLength];
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			for (int v = 0; v < 16; ++v) {
				resultVertex[v][i][j] = 0.;
			}
		}
	}
#endif

	
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		//First we construct the propagator
		matrix_t S(4*diracVectorLength,4*diracVectorLength);
		for (int c = 0; c < diracVectorLength; ++c) {
			for (int d = 0; d < diracVectorLength; ++d) {
				for (int nu = 0; nu < 4; ++nu) {
					for (int rho = 0; rho < 4; ++rho) {
						S(c*4+nu,d*4+rho) = 0.;
						for (int mu = 0; mu < 4; ++mu) {
							S(c*4+nu,d*4+rho) += vector_dot(source[c*4+nu][site][mu],inverse_source[d*4+rho][site][mu]);
						}
					}
				}
			}
		}

		for (int v = 0; v < 16; ++v) {
			matrix_t result = gamma.gamma5()*htrans(S)*gamma.gamma5()*gamma.gammaChromaMatrices(v)*S;

			for (int i = 0; i < 4*diracVectorLength; ++i) {
				for (int j = 0; j < 4*diracVectorLength; ++j) {
#ifdef MULTITHREADING
					resultVertex[v][i][j][omp_get_thread_num()] += result(i,j);
#endif
#ifndef MULTITHREADING
					resultVertex[v][i][j] += result(i,j);
#endif	
				}
			}
		}
	}
	
	//Now we collect all the results
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			//... from the threads
#ifdef MULTITHREADING
			for (int v = 0; v < 16; ++v) {
				vertex[v][i][j] = 0;
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					vertex[v][i][j] += resultVertex[v][i][j][thread];
				}
			}
#endif
#ifndef MULTITHREADING
			for (int v = 0; v < 16; ++v) {
				vertex[v][i][j] = resultVertex[v][i][j];
			}
#endif
		}
	}
	
	//... from the MPI processes
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			//Correct normalization
			long_real_t factor = 2.*diracOperator->getKappa();
			for (int v = 0; v < 16; ++v) {
				reduceAllSum(vertex[v][i][j]);
				vertex[v][i][j] = factor*factor*vertex[v][i][j]/static_cast<long_real_t>(Layout::globalVolume);
			}
		}
	}

	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("propagator");
		for (int i = 0; i < 4*diracVectorLength; ++i) {
			output->push("propagator");
			for (int j = 0; j < 4*diracVectorLength; ++j) {
				std::cout << "NPRVertex::Propagator S(" << i << "," << j << "): " << propagator[i][j] << std::endl;
				
				output->push("propagator");
				output->write("propagator", real(propagator[i][j]));
				output->write("propagator", imag(propagator[i][j]));
				output->pop("propagator");
			}
			output->pop("propagator");
		}
		output->pop("propagator");

		output->push("vertex");
		for (int v = 0; v < 16; ++v) {
			output->push("vertex");
			for (int i = 0; i < 4*diracVectorLength; ++i) {
				output->push("vertex");
				for (int j = 0; j < 4*diracVectorLength; ++j) {
					std::cout << "NPRVertex::Vertex Gamma_" << v << " V(" << i << "," << j << "): " << vertex[v][i][j] << std::endl;
				
					output->push("vertex");
					output->write("vertex", real(vertex[v][i][j]));
					output->write("vertex", imag(vertex[v][i][j]));
					output->pop("vertex");
				}
				output->pop("vertex");
			}
			output->pop("vertex");
		}
		output->pop("vertex");
	}
	
	if (redBlockDiracOperator) delete redBlockDiracOperator;
	if (blackBlockDiracOperator) delete blackBlockDiracOperator;
	if (inverter) delete inverter;
}

void NPRVertex::registerParameters(po::options_description& desc) {
	desc.add_options()
		("NPRVertex::inverter_precision", po::value<double>()->default_value(0.000000000001), "set the precision used by the inverter")
		("NPRVertex::inverter_max_steps", po::value<unsigned int>()->default_value(5000), "set the maximum steps used by the inverter")
		("NPRVertex::momentum", po::value<std::string>()->default_value("{2,2,2,2}"), "Momentum for the measure of the vertex function (syntax: {px,py,pz,pt})")
		
		("NPRVertex::multigrid", po::value<std::string>()->default_value("false"), "Should we use the multigrid inverter? true/false")
		("NPRVertex::multigrid_basis_dimension", po::value<unsigned int>()->default_value(20), "The dimension of the basis for multigrid")
		("NPRVertex::multigrid_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for Multigrid (syntax: {bx,by,bz,bt})")

		("NPRVertex::sap_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for SAP (syntax: {bx,by,bz,bt})")
		("NPRVertex::sap_iterations", po::value<unsigned int>()->default_value(5), "The number of sap iterations")
		("NPRVertex::sap_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the inner SAP inverter")
		("NPRVertex::sap_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the inner SAP inverter")
		("NPRVertex::gmres_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the GMRES inverter used to initialize the multigrid basis")
		("NPRVertex::gmres_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the GMRES inverter used to initialize the multigrid basis")
		;
}

} /* namespace Update */
