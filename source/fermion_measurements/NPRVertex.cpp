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

namespace Update {

NPRVertex::NPRVertex() : diracOperator(0), squareDiracOperator(0), biConjugateGradient(0) { }

NPRVertex::NPRVertex(const NPRVertex& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), squareDiracOperator(0), biConjugateGradient(0) { }

NPRVertex::~NPRVertex() {
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
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

	extended_fermion_lattice_t lattice;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("level_stout_smearing_meson");
		double smearingRho = environment.configurations.get<double>("rho_stout_smearing");
		StoutSmearing stoutSmearing;
#ifdef ADJOINT
		extended_gauge_lattice_t smearedConfiguration;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, smearedConfiguration, numberLevelSmearing, smearingRho);
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::convert(lattice, smearedConfiguration);//TODO
#endif
#ifndef ADJOINT
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
#endif
		environment.setFermionBc(lattice);
	} catch (NotFoundOption& ex) {
		lattice =  environment.getFermionLattice();
		if (isOutputProcess()) std::cout << "NPRVertex::No smearing options found, proceeding without!" << std::endl;
	}

	diracOperator->setLattice(lattice);
	squareDiracOperator->setLattice(lattice);
	
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));

	extended_dirac_vector_t source;

	std::complex<long_real_t> propagator[4][4];
	std::complex<long_real_t> vertex[1][4][4];

	int inversionSteps = 0;

	std::vector<real_t> momentum = environment.configurations.get< std::vector<real_t> >("npr_vertex_momentum");
	momentum[0] = momentum[0]*2.*PI/Layout::glob_x;
	momentum[1] = momentum[1]*2.*PI/Layout::glob_y;
	momentum[2] = momentum[2]*2.*PI/Layout::glob_z;
	momentum[3] = momentum[3]*2.*PI/Layout::glob_t;
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateMomentumSource(source, momentum, alpha, c);
			biConjugateGradient->solve(squareDiracOperator, source, eta);
			diracOperator->multiply(tmp[c*4 + alpha],eta);
			AlgebraUtils::gamma5(tmp[c*4 + alpha]);
			inversionSteps += biConjugateGradient->getLastSteps();
		}
	}
	
	if (isOutputProcess()) std::cout << "NPRVertex::Vertex computed with " << inversionSteps << " inversion steps" << std::endl;
	
	//Ugly initialization to zero
#ifdef MULTITHREADING
	std::complex<long_real_t> resultPropagator[4][4][omp_get_max_threads()];
	std::complex<long_real_t> resultVertex[1][4][4][omp_get_max_threads()];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < omp_get_max_threads(); ++k) {
				resultPropagator[i][j][k] = 0.;
				for (int v = 0; v < 1; ++v) {
					resultVertex[v][i][j][k] = 0.;
				}
			}
		}
	}
#endif
#ifndef MULTITHREADING
	std::complex<long_real_t> resultPropagator[4][4];
	std::complex<long_real_t> resultVertex[1][4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			resultPropagator[i][j] = 0.;
			for (int v = 0; v < 1; ++v) {
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
						resultPropagator[alpha][mu][omp_get_thread_num()] += std::complex<real_t>(cos(phase),sin(phase))*tmp[c*4+alpha][site][mu][d];
#endif
#ifndef MULTITHREADING
						resultPropagator[alpha][mu] += std::complex<real_t>(cos(phase),sin(phase))*tmp[c*4+alpha][site][mu][d];
#endif					
					}
				}
			}
		}
	}

	//Gamma5 vertex
	matrix_t gamma5(4,4);
	for (int alpha = 0; alpha < 4; ++alpha) {
		for (int beta = 0; beta < 4; ++beta) {
			if (alpha != beta) gamma5(alpha,beta) = 0.;
			else if (alpha < 2) gamma5(alpha,beta) = 1.;
			else gamma5(alpha,beta) = -1.;
		}
	}

	for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (int d = 0; d < diracVectorLength; ++d) {
				//First we construct the propagator
				matrix_t S(4,4);
				for (int nu = 0; nu < 4; ++nu) {
					for (int rho = 0; rho < 4; ++rho) {
						S(nu,rho) = tmp[c*4+nu][site][rho][d];
					}
				}
				
				matrix_t result = gamma5*htrans(S)*gamma5*gamma5*S;
				for (int nu = 0; nu < 4; ++nu) {
					for (int rho = 0; rho < 4; ++rho) {

#ifdef MULTITHREADING
						resultVertex[0][nu][rho][omp_get_thread_num()] += result(nu,rho);
#endif
#ifndef MULTITHREADING
						resultVertex[0][nu][rho] += result(nu,rho);
#endif					
					}
				}
			}
		}
	}
	
	//Now we collect all the results
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			//... from the threads
#ifdef MULTITHREADING
			propagator[i][j] = 0;
			for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
				propagator[i][j] += resultPropagator[i][j][thread];	
			}
			for (int v = 0; v < 1; ++v) {
				vertex[v][i][j] = 0;
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					vertex[v][i][j] += resultVertex[v][i][j][thread];
				}
			}
#endif
#ifndef MULTITHREADING
			propagator[i][j] = resultPropagator[i][j];
			for (int v = 0; v < 1; ++v) {
				vertex[v][i][j] = resultVertex[v][i][j][thread];
			}
#endif
		}
	}
	
	//... from the MPI processes
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			reduceAllSum(propagator[i][j].real());
			reduceAllSum(propagator[i][j].imag());
			//Correct normalization
			long_real_t factor = 2.*diracOperator->getKappa();
			propagator[i][j] = factor*propagator[i][j]/static_cast<long_real_t>(Layout::globalVolume);
			for (int v = 0; v < 1; ++v) {
				reduceAllSum(vertex[v][i][j].real());
				reduceAllSum(vertex[v][i][j].imag());
				vertex[v][i][j] = factor*factor*vertex[v][i][j]/static_cast<long_real_t>(Layout::globalVolume);
			}
		}
	}

	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("propagator");
		for (int i = 0; i < 4; ++i) {
			output->push("propagator");
			for (int j = 0; j < 4; ++j) {
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
		for (int v = 0; v < 1; ++v) {
			output->push("vertex");
			for (int i = 0; i < 4; ++i) {
				output->push("vertex");
				for (int j = 0; j < 4; ++j) {
					std::cout << "NPRVertex::Vertex " << v << " V(" << i << "," << j << "): " << vertex[v][i][j] << std::endl;
				
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

}

} /* namespace Update */
