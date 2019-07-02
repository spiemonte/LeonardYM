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
#include "inverters/PreconditionedBiCGStab.h"
#include "dirac_operators/Propagator.h"
#include "utils/MultiThreadSummator.h"

namespace Update {

MesonCorrelator::MesonCorrelator() : diracOperator(0), inverter(0) { }

MesonCorrelator::MesonCorrelator(const MesonCorrelator& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), inverter(0) { }

MesonCorrelator::~MesonCorrelator() {
	if (diracOperator) delete diracOperator;
	if (inverter) delete inverter;
}

void MesonCorrelator::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;

	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	extended_fermion_lattice_t lattice;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("MesonCorrelator::stout_smearing_levels");
		double smearingRho = environment.configurations.get<double>("MesonCorrelator::stout_smearing_rho");
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
		unsigned int t = environment.configurations.get<unsigned int>("MesonCorrelator::t_source_origin");
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
	diracOperator->setGamma5(false);
	
	if (inverter == 0) inverter = new PreconditionedBiCGStab();
	inverter->setPrecision(environment.configurations.get<double>("MesonCorrelator::inverter_precision"));
	inverter->setMaximumSteps(environment.configurations.get<unsigned int>("MesonCorrelator::inverter_max_steps"));

	std::vector< long_real_t > pionNorm;

	extended_dirac_vector_t source;

	MultiThreadSummator<long_real_t>* pionOperator = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	MultiThreadSummator<long_real_t>* scalarOperator = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	MultiThreadSummator<long_real_t>* PSOperator = new MultiThreadSummator<long_real_t>[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		pionOperator[t].reset();
		scalarOperator[t].reset();
		PSOperator[t].reset();
	}

	int inversionSteps = 0;
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateSource(source, alpha, c);
			inverter->solve(diracOperator, source, tmp);
			Propagator::constructPropagator(diracOperator, tmp, propagator[c*4 + alpha]);

			inversionSteps += inverter->getLastSteps();
		}
	}
	
	if (isOutputProcess()) std::cout << "MesonCorrelator::Correlators computed with " << inversionSteps << " inversion steps" << std::endl;
	
	//The core of the computation of the connected correlators
	//Pion correlator first
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int d = 0; d < diracVectorLength; ++d) {
						std::complex<real_t> dottmp = conj(propagator[c*4+alpha][site][mu][d])*propagator[c*4+alpha][site][mu][d];

						pionOperator[Layout::globalIndexT(site)].add(real(dottmp));
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
								std::complex<real_t> dottmp = Gamma::gamma5_gamma(3,mu,nu)*Gamma::gamma5_gamma(3,alpha,beta)*conj(propagator[c*4+alpha][site][mu][d])*propagator[c*4+beta][site][nu][d];
								scalarOperator[Layout::globalIndexT(site)].add(real(dottmp));				
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
							dottmp -= Gamma::gamma(3,alpha,beta)*conj(propagator[c*4+alpha][site][mu][d])*propagator[c*4+beta][site][mu][d];
							dottmp += Gamma::gamma(3,alpha,beta)*conj(propagator[c*4+mu][site][alpha][d])*propagator[c*4+mu][site][beta][d];

							PSOperator[Layout::globalIndexT(site)].add(real(dottmp));					
						}
					}
				}
			}
		}
	}
	
	//Now we collect all the results
	for (int t = 0; t < Layout::glob_t; ++t) {
		pionOperator[t].computeResult();
		scalarOperator[t].computeResult();
		PSOperator[t].computeResult();
	}
	
	if (environment.measurement && isOutputProcess()) {
		//Correct normalization
		real_t factor = 4.*diracOperator->getKappa()*diracOperator->getKappa();

		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("pion_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Pion Exact Correlator at t " << t << " is " << factor*pionOperator[t].getResult() << std::endl;

			output->write("pion_exact", factor*pionOperator[t].getResult());
		}
		output->pop("pion_exact");

		output->push("scalar_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Scalar Exact Correlator at t " << t << " is " << factor*scalarOperator[t].getResult() << std::endl;

			output->write("scalar_exact", factor*scalarOperator[t].getResult());
		}
		output->pop("scalar_exact");
		
		output->push("ps_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::PS Exact Correlator at t " << t << " is " << factor*PSOperator[t].getResult() << std::endl;

			output->write("ps_exact", factor*PSOperator[t].getResult());
		}
		output->pop("ps_exact");
		
		output->push("pcac_mass");
		for (int t = 0; t < Layout::glob_t; ++t) {
			int mindex = (t == 0) ? Layout::glob_t - 1 : t - 1;
			int pindex = (t == Layout::glob_t -1) ? 0 : t + 1;
			std::cout << "MesonCorrelator::PCAC mass at t " << t << " is " << (PSOperator[mindex].getResult() - PSOperator[pindex].getResult())/(4.*pionOperator[t].getResult()) << std::endl;
			
			output->write("pcac_mass", (PSOperator[mindex].getResult() - PSOperator[pindex].getResult())/(4.*pionOperator[t].getResult()));
		}
		output->pop("pcac_mass");
	}	
	
	bool compute_disconnected;
	try {
		compute_disconnected = environment.configurations.get<bool>("MesonCorrelator::compute_disconnected_contributions");
	} catch (NotFoundOption& ex) {
		compute_disconnected = false;
	}

	if (compute_disconnected) {
		MultiThreadSummator< std::complex<long_real_t> >* etaDisconnectedOperator = new MultiThreadSummator< std::complex<long_real_t> >[Layout::glob_t];
		MultiThreadSummator< std::complex<long_real_t> >* scalarDisconnectedOperator = new MultiThreadSummator< std::complex<long_real_t> >[Layout::glob_t];
		

		unsigned int numberStochasticEstimators = environment.configurations.get<unsigned int>("MesonCorrelator::number_stochastic_estimators");

		extended_dirac_vector_t randomNoise[4], inverse[4];

		for (unsigned int step = 0; step < numberStochasticEstimators; ++step) {
			for (int t = 0; t < Layout::glob_t; ++t) {
				etaDisconnectedOperator[t].reset();
				scalarDisconnectedOperator[t].reset();
			}

			for (unsigned int mu = 0; mu < 4; ++mu) {
				this->generateRandomNoise(randomNoise[mu]);
				inverter->solve(diracOperator,randomNoise[mu], inverse[mu]);
			}

#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (unsigned int nu = 0; nu < 4; ++nu) {
						std::complex<real_t> dotr = 0.;
						for (unsigned int c = 0; c < diracVectorLength; ++c) dotr += conj(randomNoise[mu][site][nu][c])*inverse[mu][site][nu][c];
						if (nu < 2) etaDisconnectedOperator[Layout::globalIndexT(site)].add(dotr);
						else etaDisconnectedOperator[Layout::globalIndexT(site)].add(-dotr);
						scalarDisconnectedOperator[Layout::globalIndexT(site)].add(dotr);
					}
				}
			}

			for (int t = 0; t < Layout::glob_t; ++t) {
				etaDisconnectedOperator[t].computeResult();
				scalarDisconnectedOperator[t].computeResult();
			}


			if (environment.measurement && isOutputProcess()) {
				long_real_t factor = 2.*diracOperator->getKappa();

				GlobalOutput* output = GlobalOutput::getInstance();

				output->push("eta_disconnected");
				for (int t = 1; t < Layout::glob_t; ++t) {
					std::cout << "Eta disconneted real contribution t " << t << " for random source " << step << " is: " << factor*etaDisconnectedOperator[t].getResult() << std::endl;

					output->write("eta_disconnected", factor*etaDisconnectedOperator[t].getResult()/static_cast<long_real_t>(Layout::glob_spatial_volume));
				}
				output->pop("eta_disconnected");

				output->push("scalar_disconnected");
				for (int t = 1; t < Layout::glob_t; ++t) {
					std::cout << "Scalar disconneted real contribution t " << t << " for random source " << step << " is: " << factor*scalarDisconnectedOperator[t].getResult() << std::endl;

					output->write("scalar_disconnected", factor*scalarDisconnectedOperator[t].getResult()/static_cast<long_real_t>(Layout::glob_spatial_volume));
				}
				output->pop("scalar_disconnected");
			}
		}

		delete[] etaDisconnectedOperator;
		delete[] scalarDisconnectedOperator;
	}

	delete[] pionOperator;
	delete[] scalarOperator;
	delete[] PSOperator;
}

void MesonCorrelator::registerParameters(po::options_description& desc) {
	static bool single = true;
	if (single) desc.add_options()
		("MesonCorrelator::inverter_precision", po::value<real_t>()->default_value(0.00000000001), "set the inverter precision")
		("MesonCorrelator::inverter_max_steps", po::value<unsigned int>()->default_value(10000), "maximum number of inverter steps")
		("MesonCorrelator::t_source_origin", po::value<unsigned int>()->default_value(0), "T origin for the wall source")
		("MesonCorrelator::stout_smearing_rho", po::value<real_t>(), "set the stout smearing parameter")
		("MesonCorrelator::stout_smearing_levels", po::value<unsigned int>(), "levels of stout smearing")
		("MesonCorrelator::number_stochastic_estimators", po::value<unsigned int>()->default_value(50), "number of stochastic estimators")
		("MesonCorrelator::compute_disconnected_contributions", po::value<bool>()->default_value(false), "Compute the disconnected contributions of the eta'?")
		;
	single = false;
}

} /* namespace Update */
