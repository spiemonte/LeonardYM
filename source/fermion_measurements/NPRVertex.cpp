/*
 * MesonCorrelator.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "NPRVertex.h"
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

NPRVertex::NPRVertex() : LatticeSweep(), MultiGridStochasticEstimator(), diracOperator(0), squareDiracOperator(0), inverter(0), gamma() { }

NPRVertex::NPRVertex(const NPRVertex& toCopy) : LatticeSweep(toCopy), MultiGridStochasticEstimator(toCopy), diracOperator(0), squareDiracOperator(0), inverter(0), gamma() { }

NPRVertex::~NPRVertex() {
	if (diracOperator) delete diracOperator;
}

void NPRVertex::execute(environment_t& environment) {
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
		GMRESR* gmresr = new GMRESR();
		inverter = gmresr;
		

		if (isOutputProcess()) std::cout << "NPRVertex::Without using multigrid ..." << std::endl;
	}

	inverter->setPrecision(environment.configurations.get<double>("NPRVertex::inverter_precision"));
	inverter->setMaximumSteps(environment.configurations.get<unsigned int>("NPRVertex::inverter_max_steps"));
	
	extended_dirac_vector_t source;

	std::complex<long_real_t> propagator[4*diracVectorLength][4*diracVectorLength];
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
			this->generateMomentumSource(source, momentum, alpha, c);
			//this->generateSource(source, alpha, c);
			inverter->solve(diracOperator, source, tmp[c*4 + alpha]);
			inversionSteps += inverter->getLastSteps();
			
			/*diracOperator->multiply(eta,source);
			biConjugateGradient->solve(squareDiracOperator, eta, tmp[c*4 + alpha]);
			AlgebraUtils::gamma5(tmp[c*4 + alpha]);
			inversionSteps += biConjugateGradient->getLastSteps();*/

			
		}
	}
	
	if (isOutputProcess()) std::cout << "NPRVertex::Vertex computed with " << inversionSteps << " inversion steps" << std::endl;
	
	//Ugly initialization to zero
#ifdef MULTITHREADING
	std::complex<long_real_t> resultPropagator[4*diracVectorLength][4*diracVectorLength][omp_get_max_threads()];
	std::complex<long_real_t> resultVertex[16][4*diracVectorLength][4*diracVectorLength][omp_get_max_threads()];
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			for (int k = 0; k < omp_get_max_threads(); ++k) {
				resultPropagator[i][j][k] = 0.;
				for (int v = 0; v < 16; ++v) {
					resultVertex[v][i][j][k] = 0.;
				}
			}
		}
	}
#endif
#ifndef MULTITHREADING
	std::complex<long_real_t> resultPropagator[4*diracVectorLength][4*diracVectorLength];
	std::complex<long_real_t> resultVertex[16][4*diracVectorLength][4*diracVectorLength];
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			resultPropagator[i][j] = 0.;
			for (int v = 0; v < 16; ++v) {
				resultVertex[v][i][j] = 0.;
			}
		}
	}
#endif
	
	//The core of the computation of the npr vertex
	//Quark propagator
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				real_t phase = - (Layout::globalIndexX(site)*momentum[0] + Layout::globalIndexY(site)*momentum[1] + Layout::globalIndexZ(site)*momentum[2] + Layout::globalIndexT(site)*momentum[3]);
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int d = 0; d < diracVectorLength; ++d) {
#ifdef MULTITHREADING
						resultPropagator[c*4 + alpha][d*4 + mu][omp_get_thread_num()] += std::complex<real_t>(cos(phase),sin(phase))*tmp[c*4+alpha][site][mu][d];
#endif
#ifndef MULTITHREADING
						resultPropagator[c*4 + alpha][d*4 + mu] += std::complex<real_t>(cos(phase),sin(phase))*tmp[c*4+alpha][site][mu][d];
#endif					
					}
				}
			}
		}
	}

	
#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		//First we construct the propagator
		matrix_t S(4*diracVectorLength,4*diracVectorLength);
		for (int c = 0; c < diracVectorLength; ++c) {
			for (int d = 0; d < diracVectorLength; ++d) {
				for (int nu = 0; nu < 4; ++nu) {
					for (int rho = 0; rho < 4; ++rho) {
						S(c*4+nu,d*4+rho) = tmp[c*4+nu][site][rho][d];
					}
				}
			}
		}

		for (int v = 0; v < 16; ++v) {
			matrix_t result = gamma.gamma5()*htrans(S)*gamma.gamma5()*gamma.gammaChromaMatrices(v)*S;
			for (int c = 0; c < diracVectorLength; ++c) {
				for (int d = 0; d < diracVectorLength; ++d) {
					for (int nu = 0; nu < 4; ++nu) {
						for (int rho = 0; rho < 4; ++rho) {
							S(c*4+nu,d*4+rho) = tmp[c*4+nu][site][rho][d];
						}
					}
				}
			}

			for (int c = 0; c < diracVectorLength; ++c) {
				for (int d = 0; d < diracVectorLength; ++d) {
					for (int nu = 0; nu < 4; ++nu) {
						for (int rho = 0; rho < 4; ++rho) {

#ifdef MULTITHREADING
							resultVertex[v][c*4+nu][d*4+rho][omp_get_thread_num()] += result(nu,rho);
#endif
#ifndef MULTITHREADING
							resultVertex[v][c*4+nu][d*4+rho] += result(nu,rho);
#endif					
						}
					}
				}
			}
		}
	}
	
	//Now we collect all the results
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			//... from the threads
#ifdef MULTITHREADING
			propagator[i][j] = 0;
			for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
				propagator[i][j] += resultPropagator[i][j][thread];	
			}
			for (int v = 0; v < 16; ++v) {
				vertex[v][i][j] = 0;
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					vertex[v][i][j] += resultVertex[v][i][j][thread];
				}
			}
#endif
#ifndef MULTITHREADING
			propagator[i][j] = resultPropagator[i][j];
			for (int v = 0; v < 16; ++v) {
				vertex[v][i][j] = resultVertex[v][i][j];
			}
#endif
		}
	}
	
	//... from the MPI processes
	for (int i = 0; i < 4*diracVectorLength; ++i) {
		for (int j = 0; j < 4*diracVectorLength; ++j) {
			reduceAllSum(propagator[i][j].real());
			reduceAllSum(propagator[i][j].imag());
			//Correct normalization
			long_real_t factor = 2.*diracOperator->getKappa();
			propagator[i][j] = factor*propagator[i][j]/static_cast<long_real_t>(Layout::globalVolume);
			for (int v = 0; v < 16; ++v) {
				reduceAllSum(vertex[v][i][j].real());
				reduceAllSum(vertex[v][i][j].imag());
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
	
	//Here we measure the disconnected contribution
	if (environment.configurations.get<std::string>("NPRVertex::measure_disconnected") == "true") {
		unsigned int numberRandomSources = environment.configurations.get<unsigned int>("NPRVertex::number_stochastic_estimators");
		extended_dirac_vector_t* randomSources = new extended_dirac_vector_t[numberRandomSources];
		extended_dirac_vector_t* inverseRandomSources = new extended_dirac_vector_t[numberRandomSources];

		MultiGridSolver* mgs = dynamic_cast<MultiGridSolver*>(inverter);
		if (mgs && environment.configurations.get<std::string>("NPRVertex::use_multigrid_diluition") == "true") {
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
			if (isOutputProcess()) std::cout << "NPRVertex::Disconnected contributions measured without multigrid" << std::endl;
			for (unsigned int i = 0; i < numberRandomSources; ++i) {
				this->generateRandomNoise(randomSources[i]);
				inverter->solve(diracOperator, randomSources[i], inverseRandomSources[i]);
			}
		}

		HoppingOperator H(diracOperator);
		GammaOperators gammaOperators;
		extended_dirac_vector_t tmp;
		/*extended_dirac_vector_t tmp, tmp2, tmp3, tmp4;
		H.apply(tmp, inverseRandomSources[0]);
		H.apply(tmp2, tmp);
		diracOperator->multiply(tmp, inverseRandomSources[0]);
		diracOperator->multiply(tmp3, tmp);
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmp4[site][mu] = tmp3[site][mu] - 2*tmp[site][mu] + inverseRandomSources[0][site][mu];
			}
		}
		std::cout << "*************** Giusto per AAA+: " << AlgebraUtils::differenceNorm(tmp4, tmp2) << std::endl;
		std::cout << "*************** Giusto per AAA-: " << AlgebraUtils::differenceNorm(tmp, randomSources[0]) << std::endl;*/
		

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

		extended_dirac_vector_t source, hopping_source, previous;
				
		for (unsigned int i = 0; i < 3000; ++i) {
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
	}

	
	if (redBlockDiracOperator) delete redBlockDiracOperator;
	if (blackBlockDiracOperator) delete blackBlockDiracOperator;
	if (inverter) delete inverter;
}

void NPRVertex::registerParameters(po::options_description& desc) {
	desc.add_options()
		("NPRVertex::inverter_precision", po::value<double>()->default_value(0.0000000001), "set the precision used by the inverter")
		("NPRVertex::inverter_max_steps", po::value<unsigned int>()->default_value(5000), "set the maximum steps used by the inverter")
		("NPRVertex::momentum", po::value<std::string>()->default_value("{2,2,2,2}"), "Momentum for the measure of the vertex function (syntax: {px,py,pz,pt})")
		
		("NPRVertex::multigrid", po::value<std::string>()->default_value("true"), "Should we use the multigrid inverter? true/false")
		("NPRVertex::multigrid_basis_dimension", po::value<unsigned int>()->default_value(20), "The dimension of the basis for multigrid")
		("NPRVertex::multigrid_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for Multigrid (syntax: {bx,by,bz,bt})")

		("NPRVertex::sap_block_size", po::value<std::string>()->default_value("{4,4,4,4}"), "Block size for SAP (syntax: {bx,by,bz,bt})")
		("NPRVertex::sap_iterations", po::value<unsigned int>()->default_value(5), "The number of sap iterations")
		("NPRVertex::sap_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the inner SAP inverter")
		("NPRVertex::sap_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the inner SAP inverter")
		("NPRVertex::gmres_inverter_precision", po::value<double>()->default_value(0.00000000001), "The precision of the GMRES inverter used to initialize the multigrid basis")
		("NPRVertex::gmres_inverter_max_steps", po::value<unsigned int>()->default_value(100), "The maximum number of steps for the GMRES inverter used to initialize the multigrid basis")
		
		("NPRVertex::measure_disconnected", po::value<std::string>()->default_value("true"), "Should we measure the disconnected contributions (default: true)")
		("NPRVertex::number_stochastic_estimators", po::value<unsigned int>()->default_value(13), "Number of stochastic estimators for the disconnected part")
		("NPRVertex::use_multigrid_diluition", po::value<std::string>()->default_value("true"), "Should we use the multigrid diluition for the measure of disconnected contribution? true/false")
		;
}

} /* namespace Update */
