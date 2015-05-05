/*
 * MesonCorrelator.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "MesonCorrelator.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/StoutSmearing.h"
#include "utils/Gamma.h"
#ifndef PI
#define PI 3.141592653589793238462643
#endif

namespace Update {

MesonCorrelator::MesonCorrelator() : diracOperator(0), squareDiracOperator(0), biConjugateGradient(0) { }

MesonCorrelator::MesonCorrelator(const MesonCorrelator& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), squareDiracOperator(0), biConjugateGradient(0) { }

MesonCorrelator::~MesonCorrelator() {
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
}

void MesonCorrelator::execute(environment_t& environment) {
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
		if (isOutputProcess()) std::cout << "MesonCorrelator::No smearing options found, proceeding without!" << std::endl;
	}

	try {
		unsigned int t = environment.configurations.get<unsigned int>("t_source_origin");
		extended_fermion_lattice_t swaplinkconfig;
		typedef extended_fermion_lattice_t LT;
		if (t != 0) {
			for (unsigned int n = 0; n < t; ++n) {
				//We do a swap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = lattice[site][mu];
				}
				swaplinkconfig.updateHalo();
				//We wrap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) lattice[site][mu] = swaplinkconfig[LT::sup(site,3)][mu];
				}
				lattice.updateHalo();
			}
		}
	} catch (NotFoundOption& ex) {

	}

	diracOperator->setLattice(lattice);
	squareDiracOperator->setLattice(lattice);
	
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));

	std::vector< long_real_t > pionNorm;

	extended_dirac_vector_t source;

	long_real_t pionOperator[Layout::glob_t];
	long_real_t scalarOperator[Layout::glob_t];
	long_real_t PSOperator[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		pionOperator[t] = 0.;
		scalarOperator[t] = 0.;
		PSOperator[t] = 0.;
	}

	int inversionSteps = 0;
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateSource(source, alpha, c);
			biConjugateGradient->solve(squareDiracOperator, source, eta);
			diracOperator->multiply(tmp[c*4 + alpha],eta);
			inversionSteps += biConjugateGradient->getLastSteps();
		}
	}
	
	if (isOutputProcess()) std::cout << "MesonCorrelator::Correlators computed with " << inversionSteps << " inversion steps" << std::endl;
	
	//Ugly initialization to zero
#ifdef MULTITHREADING
	long_real_t resultPion[Layout::glob_t][omp_get_max_threads()];
	long_real_t resultScalar[Layout::glob_t][omp_get_max_threads()];
	long_real_t resultPS[Layout::glob_t][omp_get_max_threads()];
	for (int i = 0; i < Layout::glob_t; ++i) {
		for (int j = 0; j < omp_get_max_threads(); ++j) {
			resultPion[i][j] = 0.;
			resultScalar[i][j] = 0.;
			resultPS[i][j] = 0.;
		}
	}
#endif
#ifndef MULTITHREADING
	long_real_t resultPion[Layout::glob_t];
	long_real_t resultScalar[Layout::glob_t];
	long_real_t resultPS[Layout::glob_t];
	for (int i = 0; i < Layout::glob_t; ++i) {
		resultPion[i] = 0.;
		resultScalar[i] = 0.;
		resultPS[i] = 0.;
	}
#endif
	
	//The core of the computation of the connected correlators
	//Pion correlator first
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int d = 0; d < diracVectorLength; ++d) {
						std::complex<real_t> dottmp = conj(tmp[c*4+alpha][site][mu][d])*tmp[c*4+alpha][site][mu][d];
#ifdef MULTITHREADING
						resultPion[Layout::globalIndexT(site)][omp_get_thread_num()] += real(dottmp);
#endif
#ifndef MULTITHREADING
						resultPion[Layout::globalIndexT(site)] += real(dottmp);
#endif					
					}
				}
			}
		}
	}
	
	//Scalar connected correlator then
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (unsigned int beta = 0; beta < 4; ++beta) {
			for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
				for (int site = 0; site < Layout::localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						for (unsigned int nu = 0; nu < 4; ++nu) {
							for (int d = 0; d < diracVectorLength; ++d) {
								std::complex<real_t> dottmp = Gamma::gamma5_gamma(3,mu,nu)*Gamma::gamma5_gamma(3,alpha,beta)*conj(tmp[c*4+alpha][site][mu][d])*tmp[c*4+beta][site][nu][d];
#ifdef MULTITHREADING
								resultScalar[Layout::globalIndexT(site)][omp_get_thread_num()] += real(dottmp);
#endif
#ifndef MULTITHREADING
								resultScalar[Layout::globalIndexT(site)] += real(dottmp);
#endif						
							}
						}
					}
				}
			}
		}
	}
	
	//PseudoScalar-PseudoVector connected correlator then
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
	for (unsigned int beta = 0; beta < 4; ++beta) {
		for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					
						for (int d = 0; d < diracVectorLength; ++d) {
							std::complex<real_t> dottmp(0.,0.);
							dottmp -= Gamma::gamma(3,alpha,beta)*conj(tmp[c*4+alpha][site][mu][d])*tmp[c*4+beta][site][mu][d];
							dottmp += Gamma::gamma(3,alpha,beta)*conj(tmp[c*4+mu][site][alpha][d])*tmp[c*4+mu][site][beta][d];
#ifdef MULTITHREADING
							resultPS[Layout::globalIndexT(site)][omp_get_thread_num()] += real(dottmp);
#endif
#ifndef MULTITHREADING
							resultPS[Layout::globalIndexT(site)] += real(dottmp);
#endif						
						}
					}
				}
			}
		}
	}
	
	//Now we collect all the results
	for (int t = 0; t < Layout::glob_t; ++t) {
		//... from the threads
#ifdef MULTITHREADING
		for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
			pionOperator[t] += resultPion[t][thread];
			scalarOperator[t] += resultScalar[t][thread];
			PSOperator[t] += resultPS[t][thread];
		}
#endif
#ifndef MULTITHREADING
		pionOperator[t] += resultPion[t];
		scalarOperator[t] += resultScalar[t];
		PSOperator[t] += resultPS[t];
#endif
	}
	
	//... from the MPI processes
	for (int t = 0; t < Layout::glob_t; ++t) {
		reduceAllSum(pionOperator[t]);
		reduceAllSum(scalarOperator[t]);
		reduceAllSum(PSOperator[t]);
		//Correct normalization
		real_t factor = 4.*diracOperator->getKappa()*diracOperator->getKappa();
		pionOperator[t] = factor*pionOperator[t];
		scalarOperator[t] = factor*scalarOperator[t];
		PSOperator[t] = factor*PSOperator[t];
	}

	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("pion_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Pion Exact Correlator at t " << t << " is " << pionOperator[t] << std::endl;

			output->write("pion_exact", pionOperator[t]);
		}
		output->pop("pion_exact");

		output->push("scalar_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Scalar Exact Correlator at t " << t << " is " << scalarOperator[t] << std::endl;

			output->write("scalar_exact", scalarOperator[t]);
		}
		output->pop("scalar_exact");
		
		output->push("ps_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::PS Exact Correlator at t " << t << " is " << PSOperator[t] << std::endl;

			output->write("ps_exact", PSOperator[t]);
		}
		output->pop("ps_exact");
		
		output->push("pcac_mass");
		for (int t = 0; t < Layout::glob_t; ++t) {
			int mindex = (t == 0) ? Layout::glob_t - 1 : t - 1;
			int pindex = (t == Layout::glob_t -1) ? 0 : t + 1;
			std::cout << "MesonCorrelator::PCAC mass at t " << t << " is " << (PSOperator[mindex] - PSOperator[pindex])/(4.*pionOperator[t]) << std::endl;
			
			output->write("pcac_mass", (PSOperator[mindex] - PSOperator[pindex])/(4.*pionOperator[t]));
		}
		output->pop("pcac_mass");
	}


/*
	for (int t = 0; t < Layout::glob_t; ++t) pionOperator2[t] = 0.;

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoiseD, 0);
		for (unsigned int mu = 0; mu < 4; ++mu) {
			biConjugateGradient->solve(squareDiracOperator,randomNoiseD[mu], tmpb);
			diracOperator->multiply(tmpD[mu],tmpb);
		}

#ifdef MULTITHREADING
		long_real_t result[Layout::glob_t][omp_get_max_threads()];
		for (int i = 0; i < Layout::glob_t; ++i) {
			for (int j = 0; j < omp_get_max_threads(); ++j) {
				result[i][j] = 0.;
			}
		}
#endif
#ifndef MULTITHREADING
		long_real_t result[Layout::glob_t];
		for (int i = 0; i < Layout::glob_t; ++i) {
			result[i] = 0.;
		}
#endif

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
					result[Layout::globalIndexT(site)][omp_get_thread_num()] += real(vector_dot(tmpD[nu][site][mu],tmpD[nu][site][mu]));
#endif
#ifndef MULTITHREADING
					result[Layout::globalIndexT(site)] += real(vector_dot(tmpD[nu][site][mu],tmpD[nu][site][mu]));
#endif					
				}
			}
		}

		for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
			pionOperator2[t] = 0;
			for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
				pionOperator2[t] += result[t][thread];
			}
#endif
#ifndef MULTITHREADING
			pionOperator2[t] = result[t];
#endif
		}
	}		

	for (int t = 0; t < Layout::glob_t; ++t) {
		reduceAllSum(pionOperator2[t]);
	}

	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("pion_stoch");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "Pion Stochastic Correlator at t " << t << " is " << pionOperator2[t] << std::endl;

			output->write("pion_stoch", pionOperator2[t]);
		}
		output->pop("pion_stoch");
	}

*/



	


	
	
	/*bool calculate_disconnected;
	try {
		calculate_disconnected = environment.configurations.get<bool>("calculate_disconnected_contributions");
	} catch (NotFoundOption& ex) {
		calculate_disconnected = false;
	}

	if (calculate_disconnected) {
		for (unsigned int t = 0; t < Layout::glob_t; ++t) {
			for (unsigned int step = 0; step < max_step; ++step) {
				this->generateRandomNoise(randomNoiseD, t);
				for (unsigned int mu = 0; mu < 4; ++mu) {
					biConjugateGradient->solve(diracOperator,randomNoiseD[mu], tmpD[mu]);
				}

				//Local eta operator
				long_real_t result_re = 0.;

				//only for this t
				if (t/Layout::loc_t == Layout::rankT()) {
					int local_t = (t % Layout::loc_t);
					//calculate the result
#ifdef MULTITHREADING
#pragma omp parallel for reduction(+:result_re)
#endif
					for (int x = 0; x < Layout::loc_x; ++x) {
						for (int y = 0; y < Layout::loc_y; ++y) {
							for (int z = 0; z < Layout::loc_z; ++z) {
								int site = Layout::index(x,y,z,local_t);
								for (unsigned int mu = 0; mu < 4; ++mu) {
									result_re += real(vector_dot(randomNoiseD[mu][site][mu],tmpD[mu][site][mu]));
								}
							}
						}
					}
				}
#ifdef ENABLE_MPI
				::MpiExt::reduceAllSum(result_re);
#endif
				etaOperatorRe[t].push_back(result_re);
				//
				//Layout::tIndex(site);
				//Layout::toGlobalCoordinateT(t)
				//reduceAllSumArray
			}
		}

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("eta_disconnected");
			for (int t = 1; t < Layout::glob_t; ++t) {
				std::cout << "Eta disconneted real contribution diluited at t " << t << " is " << this->mean(etaOperatorRe[t]) << " +/- " << this->standardDeviation(etaOperatorRe[t]) << std::endl;

				output->write("eta_disconnected", this->mean(etaOperatorRe[t]));
			}
			output->pop("eta_disconnected");
		}
	}*/
	//delete[] pionOperator;
	//delete[] etaOperatorRe;
	//delete[] etaOperatorIm;
}

} /* namespace Update */
