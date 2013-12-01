#include "GluinoGlue.h"
#include "LieGenerators.h"
#include "Gamma.h"
#include "StoutSmearing.h"
#include "GlobalOutput.h"
#include "AlgebraUtils.h"

namespace Update {

GluinoGlue::GluinoGlue() : LatticeSweep(), StochasticEstimator(), diracOperator(0), biConjugateGradient(0) { }

GluinoGlue::GluinoGlue(const GluinoGlue& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), biConjugateGradient(0) { }

GluinoGlue::~GluinoGlue() {
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
}

void GluinoGlue::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;
	LieGenerator<GaugeGroup> tau;

	unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");

	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	extended_fermion_lattice_t lattice;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("level_stout_smearing_gluinoglue");
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
		switchAntiperiodicBc(lattice);
	} catch (NotFoundOption& ex) {
		lattice = environment.getFermionLattice();
		if (isOutputProcess()) std::cout << "GluinoGlue::No smearing options found, proceeding without!" << std::endl;
	}

	diracOperator->setLattice(lattice);
	diracOperator->setGamma5(false);

	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("inverter_max_steps"));
	
	long_real_t gluinoGlueCorrelator[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		gluinoGlueCorrelator[t] = 0.;
	}

	for (int t = 0; t < Layout::glob_t; ++t) {
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			//First we generate the source and eta
#pragma omp parallel for
			for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						eta[site][mu][c] = 0;
						source[site][mu][c] = 0.;
						if (Layout::globalIndexT(site) == 0) {
							for (unsigned int i = 0; i < 3; ++i) {
								for (unsigned int j = 0; j < 3; ++j) {
									if (Sigma::sigma(i,j,alpha,mu) != 0.) eta[site][mu][c] += Sigma::sigma(i,j,alpha,mu)*trace(cloverPlaquette(environment.gaugeLinkConfiguration,site,i,j)*tau.get(c));
								}
							}
						}
						if (Layout::globalIndexT(site) == t) {
							for (unsigned int i = 0; i < 3; ++i) {
								for (unsigned int j = 0; j < 3; ++j) {
									if (Sigma::sigma(i,j,alpha,mu) != 0.) source[site][mu][c] += Sigma::sigma(i,j,alpha,mu)*trace(cloverPlaquette(environment.gaugeLinkConfiguration,site,i,j)*tau.get(c));
								}
							}
						}
					}
				}
			}

			biConjugateGradient->solve(diracOperator, source, psi);

			long_real_t result = 0.;
#pragma omp parallel for reduction(+:result)
			for (int site = 0; site < Layout::localsize; ++site) {
				if (Layout::globalIndexT(site) == 0) {
					for (unsigned int beta = 0; beta < 4; ++beta) {
						for (int a = 0; a < diracVectorLength; ++a) {
							result += real(eta[site][beta][a]*psi[site][beta][a]);
						}
					}
				}
			}
			gluinoGlueCorrelator[t] = result;
/*
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
				std::complex<real_t> tmp(0.,0.);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					for (int a = 0; a < diracVectorLength; ++a) {
						tmp += eta[site][beta][a]*psi[site][beta][a];
					}
				}
#ifdef MULTITHREADING
				result[Layout::globalIndexT(site)][omp_get_thread_num()] += real(tmp);
#endif
#ifndef MULTITHREADING
				result[Layout::globalIndexT(site)] += real(tmp);
#endif
			}

			//We collect the results
			for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					gluinoGlueCorrelator[t] += result[t][thread];
				}
#endif
#ifndef MULTITHREADING
				gluinoGlueCorrelator[t] += result[t];
#endif
			}*/
		}
	}
	/*for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		//First we generate the source and eta
#pragma omp parallel for
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					eta[alpha][site][mu][c] = 0;
					for (unsigned int i = 0; i < 3; ++i) {
						for (unsigned int j = 0; j < 3; ++j) {
							eta[alpha][site][mu][c] += Sigma::sigma(i,j,alpha,mu)*trace(cloverPlaquette(environment.gaugeLinkConfiguration,site,i,j)*tau.get(c));
						}
					}
					if (Layout::globalIndexT(site) == 0) {
						psi[alpha][site][mu][c] = eta[alpha][site][mu][c];
					}
					else {
						psi[alpha][site][mu][c] = 0.;
					}
				}
			}
		}
	}

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise);
		biConjugateGradient->solve(diracOperator, randomNoise, rho);
		
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			std::complex<real_t> projection = static_cast< std::complex<real_t> >(AlgebraUtils::dot(randomNoise,psi[alpha]));
			
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
				std::complex<real_t> tmp(0.,0.);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					for (int a = 0; a < diracVectorLength; ++a) {
						tmp += conj(eta[alpha][site][beta][a])*rho[site][beta][a]*projection;
					}
				}
#ifdef MULTITHREADING
				result[Layout::globalIndexT(site)][omp_get_thread_num()] += real(tmp);
				if (site == 0) std::cout << "Giusto per: " << imag(tmp) << std::endl;
#endif
#ifndef MULTITHREADING
				result[Layout::globalIndexT(site)] += real(tmp);
#endif
			}

			//We collect the results
			for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					gluinoGlueCorrelator[t] += result[t][thread];
				}
#endif
#ifndef MULTITHREADING
				gluinoGlueCorrelator[t] += result[t];
#endif
			}
		}
	}*/

	for (int t = 0; t < Layout::glob_t; ++t) {
		reduceAllSum(gluinoGlueCorrelator[t]);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("gluinoglue");
		for (int t = 0; t < Layout::glob_t; ++t) {
			std::cout << "GluinoGlue::Correlator at t " << t << " is " << -gluinoGlueCorrelator[t]/max_step << std::endl;

			output->write("gluinoglue", -gluinoGlueCorrelator[t]);
		}
		output->pop("gluinoglue");
	}
	
}

GaugeGroup GluinoGlue::cloverPlaquette(const extended_gauge_lattice_t& lattice, int site, int mu, int nu) {
	typedef extended_dirac_vector_t LT;
	
	GaugeGroup result;
	set_to_zero(result);
	if (mu != nu) {
		result += lattice[site][mu]*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
		result += htrans(lattice[LT::sdn(site,nu)][nu])*lattice[LT::sdn(site,nu)][mu]*lattice[LT::sdn(LT::sup(site,mu),nu)][nu]*htrans(lattice[site][mu]);
		result += lattice[site][nu]*htrans(lattice[LT::sdn(LT::sup(site,nu),mu)][mu])*htrans(lattice[LT::sdn(site,mu)][nu])*lattice[LT::sdn(site,mu)][mu];
		result += htrans(lattice[LT::sdn(site,mu)][mu])*htrans(lattice[LT::sdn(LT::sdn(site,mu),nu)][nu])*lattice[LT::sdn(LT::sdn(site,mu),nu)][mu]*lattice[LT::sdn(site,nu)][nu];
	}
	GaugeGroup final = std::complex<real_t>(0,+1./8.)*(result-htrans(result));
	return final;
}

}

