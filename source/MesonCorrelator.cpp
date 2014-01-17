/*
 * MesonCorrelator.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "MesonCorrelator.h"
#include "BiConjugateGradient.h"
#include "GlobalOutput.h"
#include "AlgebraUtils.h"
#include "StoutSmearing.h"
#ifndef PI
#define PI 3.141592653589793238462643
#endif

namespace Update {

MesonCorrelator::MesonCorrelator() : diracOperator(0), biConjugateGradient(0) { }

MesonCorrelator::MesonCorrelator(const MesonCorrelator& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), biConjugateGradient(0) { }

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
		switchAntiperiodicBc(lattice);
	} catch (NotFoundOption& ex) {
		lattice =  environment.getFermionLattice();
		if (isOutputProcess()) std::cout << "MesonCorrelator::No smearing options found, proceeding without!" << std::endl;
	}

	diracOperator->setLattice(lattice);
	diracOperator->setGamma5(false);
	
	if (biConjugateGradient == 0) biConjugateGradient = new BiConjugateGradient();
	biConjugateGradient->setPrecision(environment.configurations.get<double>("generic_inverter_precision"));
	biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));

	std::vector< long_real_t > pionNorm;
	//std::vector< long_real_t >* pionOperator = new std::vector< long_real_t >[Layout::pgrid_t*Layout::loc_t];
	std::vector< long_real_t >* etaOperatorRe = new std::vector< long_real_t >[Layout::pgrid_t*Layout::loc_t];
	std::vector< long_real_t >* etaOperatorIm = new std::vector< long_real_t >[Layout::pgrid_t*Layout::loc_t];

	extended_dirac_vector_t source;
	extended_dirac_vector_t tmp, tmpb;

	real_t momenta[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		momenta[t] = 2.*PI*(1.-0.5*Layout::glob_t+static_cast<real_t>(t))/static_cast<real_t>(Layout::glob_t);
	}

	real_t momenta_x1 = 2.*PI*(static_cast<real_t>(4))/static_cast<real_t>(Layout::glob_t);
	real_t momenta_y1 = 2.*PI*(static_cast<real_t>(3))/static_cast<real_t>(Layout::glob_t);

	real_t momenta_x2 = 2.*PI*(static_cast<real_t>(5))/static_cast<real_t>(Layout::glob_t);
	real_t momenta_y2 = 2.*PI*(static_cast<real_t>(0))/static_cast<real_t>(Layout::glob_t);
	
	long_real_t pionOperatorZero[Layout::glob_t];
	long_real_t vectorOperatorZero[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		pionOperatorZero[t] = 0.;
		vectorOperatorZero[t] = 0.;
	}

	long_real_t pionOperatorMomentaAsymmetry1[Layout::glob_t];
	long_real_t pionOperatorMomentaAsymmetry2[Layout::glob_t];
	for (int t = 0; t < Layout::glob_t; ++t) {
		pionOperatorMomentaAsymmetry1[t] = 0.;
		pionOperatorMomentaAsymmetry2[t] = 0.;
	}

	long_real_t pionOperatorMomenta[Layout::glob_t][Layout::glob_t];
	for (int m = 0; m < Layout::glob_t; ++m) {
		for (int t = 0; t < Layout::glob_t; ++t) pionOperatorMomenta[m][t] = 0.;
	}

	double check1 = 0., check2 = 0.;

	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < numberColors; ++c) {
			this->generateSource(source, alpha, c);
			biConjugateGradient->solve(diracOperator, source, tmp);
			
#ifdef MULTITHREADING
			long_real_t resultZero[Layout::glob_t][omp_get_max_threads()];
			long_real_t resultVector[Layout::glob_t][omp_get_max_threads()];
			long_real_t resultMomentaAsymmetry1[Layout::glob_t][omp_get_max_threads()];
			long_real_t resultMomentaAsymmetry2[Layout::glob_t][omp_get_max_threads()];
			long_real_t resultMomenta[Layout::glob_t][Layout::glob_t][omp_get_max_threads()];
			for (int i = 0; i < Layout::glob_t; ++i) {
				for (int j = 0; j < omp_get_max_threads(); ++j) {
					resultZero[i][j] = 0.;
					resultVector[i][j] = 0.;
					resultMomentaAsymmetry1[i][j] = 0.;
					resultMomentaAsymmetry2[i][j] = 0.;
					for (int m = 0; m < Layout::glob_t; ++m) {
						resultMomenta[m][i][j] = 0.;
					}
				}
			}
#endif
#ifndef MULTITHREADING
			long_real_t resultZero[Layout::glob_t];
			long_real_t resultVector[Layout::glob_t];
			long_real_t resultMomentaAsymmetry1[Layout::glob_t];
			long_real_t resultMomentaAsymmetry2[Layout::glob_t];
			long_real_t resultMomenta[Layout::glob_t][Layout::glob_t];
			for (int i = 0; i < Layout::glob_t; ++i) {
				resultZero[i] = 0.;
				resultVector[i] = 0.;
				resultMomentaAsymmetry1[i] = 0.;
				resultMomentaAsymmetry2[i] = 0.;
				for (int m = 0; m < Layout::glob_t; ++m) {
					resultMomenta[m][i] = 0.;
				}
			}
#endif
			
#pragma omp parallel for
			for (int site = 0; site < Layout::localsize; ++site) {
				//Pion correlator first
				for (unsigned int mu = 0; mu < 4; ++mu) {
					std::complex<real_t> dottmp = vector_dot(tmp[site][mu],tmp[site][mu]);
#ifdef MULTITHREADING
					resultZero[Layout::globalIndexT(site)][omp_get_thread_num()] += real(dottmp);
					resultMomentaAsymmetry1[Layout::globalIndexT(site)][omp_get_thread_num()] += real(std::complex<real_t>(cos(momenta_x1*Layout::globalIndexX(site)+momenta_y1*Layout::globalIndexY(site)),sin(momenta_x1*Layout::globalIndexX(site)+momenta_y1*Layout::globalIndexY(site)))*dottmp);
					resultMomentaAsymmetry2[Layout::globalIndexT(site)][omp_get_thread_num()] += real(std::complex<real_t>(cos(momenta_x2*Layout::globalIndexX(site)+momenta_y2*Layout::globalIndexY(site)),sin(momenta_x2*Layout::globalIndexX(site)+momenta_y2*Layout::globalIndexY(site)))*dottmp);
					for (int m = 0; m < Layout::glob_t; ++m) {
						resultMomenta[m][Layout::globalIndexT(site)][omp_get_thread_num()] += real(std::complex<real_t>(cos(momenta[m]*Layout::globalIndexX(site)),sin(momenta[m]*Layout::globalIndexX(site)))*dottmp);
					}
#endif
#ifndef MULTITHREADING
					resultZero[Layout::globalIndexT(site)] += real(dottmp);
					resultMomentaAsymmetry1[Layout::globalIndexT(site)] += real(std::complex<real_t>(cos(momenta_x1*Layout::globalIndexX(site)+momenta_y1*Layout::globalIndexY(site)),sin(momenta_x1*Layout::globalIndexX(site)+momenta_y1*Layout::globalIndexY(site)))*dottmp);
					resultMomentaAsymmetry2[Layout::globalIndexT(site)] += real(std::complex<real_t>(cos(momenta_x2*Layout::globalIndexX(site)+momenta_y2*Layout::globalIndexY(site)),sin(momenta_x2*Layout::globalIndexX(site)+momenta_y2*Layout::globalIndexY(site)))*dottmp);
					for (int m = 0; m < Layout::glob_t; ++m) {
						resultMomenta[m][Layout::globalIndexT(site)] += real(std::complex<real_t>(cos(momenta[m]*Layout::globalIndexX(site)),sin(momenta[m]*Layout::globalIndexX(site)))*dottmp);
					}
#endif					
				}
				//Vector correlator then
				std::complex<real_t> dottmp = - vector_dot(tmp[site][0],tmp[site][0]) + vector_dot(tmp[site][1],tmp[site][1]) - vector_dot(tmp[site][2],tmp[site][2]) + vector_dot(tmp[site][3],tmp[site][3]);
#ifdef MULTITHREADING
				resultVector[Layout::globalIndexT(site)][omp_get_thread_num()] += real(dottmp);
#endif
#ifndef MULTITHREADING
				resultVector[Layout::globalIndexT(site)] += real(dottmp);
#endif
			}
			
			for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
				for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
					pionOperatorZero[t] += resultZero[t][thread];
					vectorOperatorZero[t] += resultVector[t][thread];
					pionOperatorMomentaAsymmetry1[t] += resultMomentaAsymmetry1[t][thread];
					pionOperatorMomentaAsymmetry2[t] += resultMomentaAsymmetry2[t][thread];
					for (int m = 0; m < Layout::glob_t; ++m) {
						pionOperatorMomenta[m][t] += resultMomenta[m][t][thread];
					}
				}
#endif
#ifndef MULTITHREADING
				pionOperatorZero[t] += resultZero[t];
				vectorOperatorZero[t] += resultVector[t];
				pionOperatorMomentaAsymmetry1[t] += resultMomentaAsymmetry1[t];
				pionOperatorMomentaAsymmetry2[t] += resultMomentaAsymmetry2[t];
				for (int m = 0; m < Layout::glob_t; ++m) {
					pionOperatorMomenta[m][t] += resultMomenta[m][t];
				}
#endif
			}

			Lattice::Site site1(0,0,0,5);
			Lattice::Site site2(3,0,0,4);

			if (Layout::localIndex[Layout::getGlobalCoordinate(site1)] != -1) {
				int site = Layout::localIndex[Layout::getGlobalCoordinate(site1)];
				for (unsigned int mu = 0; mu < 4; ++mu) {
					check1 += real(vector_dot(tmp[site][mu],tmp[site][mu]));
				}
			}

			if (Layout::localIndex[Layout::getGlobalCoordinate(site2)] != -1) {
				int site = Layout::localIndex[Layout::getGlobalCoordinate(site2)];
				for (unsigned int mu = 0; mu < 4; ++mu) {
					check2 += real(vector_dot(tmp[site][mu],tmp[site][mu]));
				}
			}
		}
	}

	reduceAllSum(check1);
	reduceAllSum(check2);
	
	for (int t = 0; t < Layout::glob_t; ++t) {
		reduceAllSum(pionOperatorZero[t]);
		reduceAllSum(vectorOperatorZero[t]);
		reduceAllSum(pionOperatorMomentaAsymmetry1[t]);
		reduceAllSum(pionOperatorMomentaAsymmetry2[t]);
		for (int m = 0; m < Layout::glob_t; ++m) {
			reduceAllSum(pionOperatorMomenta[m][t]);
		}
	}

	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("pion_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Pion Exact Correlator at t " << t << " is " << pionOperatorZero[t] << std::endl;

			output->write("pion_exact", pionOperatorZero[t]);
		}
		output->pop("pion_exact");

		output->push("vector_exact");
		for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
			std::cout << "MesonCorrelator::Vector Exact Correlator at t " << t << " is " << vectorOperatorZero[t] << std::endl;

			output->write("vector_exact", vectorOperatorZero[t]);
		}
		output->pop("vector_exact");

		output->push("pion_momenta");
		for (int m = 0; m < Layout::glob_t; ++m) {
			output->push("pion_momenta");
			for (int t = 0; t < Layout::glob_t; ++t) {
				output->write("pion_momenta", pionOperatorMomenta[m][t]);
			}
			output->pop("pion_momenta");
		}
		output->pop("pion_momenta");

		output->push("pion_momenta_asymmetry1");
		for (int t = 0; t < Layout::glob_t; ++t) {
			output->write("pion_momenta_asymmetry1", pionOperatorMomentaAsymmetry1[t]);
		}
		output->pop("pion_momenta_asymmetry1");
		output->push("pion_momenta_asymmetry2");
		for (int t = 0; t < Layout::glob_t; ++t) {
			output->write("pion_momenta_asymmetry2", pionOperatorMomentaAsymmetry2[t]);
		}
		output->pop("pion_momenta_asymmetry2");

		output->push("pion_asymmetry");
		output->write("pion_asymmetry", check1);
		output->write("pion_asymmetry", check2);
		output->pop("pion_asymmetry");
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
	delete[] etaOperatorRe;
	delete[] etaOperatorIm;
}

} /* namespace Update */
