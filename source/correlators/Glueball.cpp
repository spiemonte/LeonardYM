/*
 * Glueball.cpp
 *
 *  Created on: Oct 9, 2012
 *      Author: spiem_01
 */

#include "Glueball.h"
#include "io/GlobalOutput.h"
#include "utils/ToString.h"
#include "utils/StoutSmearing.h"
#ifndef PI
#define PI 3.14159265358979323846264338327950288419
#endif

namespace Update {

Glueball::Glueball() { }

Glueball::~Glueball() { }

void Glueball::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	extended_gauge_lattice_t lattice;
	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("level_stout_smearing_glueball");
		double smearingRho = environment.configurations.get<double>("rho_stout_smearing");
		extended_gauge_lattice_t smearedConfiguration;
		StoutSmearing stoutSmearing;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
	} catch (NotFoundOption& ex) {
		lattice =  environment.gaugeLinkConfiguration;
		if (isOutputProcess()) std::cout << "Glueball::No smearing options found, proceeding without!" << std::endl;
	}

	//First we measure the zero momentum projection
	{
		long_real_t zero_glueball[Layout::glob_t];
		long_real_t two_glueball[Layout::glob_t];
#ifdef MULTITHREADING
		long_real_t result_zero[Layout::glob_t][omp_get_max_threads()];
		long_real_t result_two[Layout::glob_t][omp_get_max_threads()];
		for (int i = 0; i < Layout::glob_t; ++i) {
			for (int j = 0; j < omp_get_max_threads(); ++j) {
				result_zero[i][j] = 0.;
				result_two[i][j] = 0.;
			}
		}
#endif
#ifndef MULTITHREADING
		long_real_t result_zero[Layout::glob_t];
		long_real_t result_two[Layout::glob_t];
		for (int i = 0; i < Layout::glob_t; ++i) {
			result_zero[i] = 0.;
			result_two[i] = 0.;
		}
#endif

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			//long_real_t result_phi_x = real(trace((lattice[LT::sup(site, 0)][1])*htrans(lattice[LT::sup(site, 1)][0])*htrans(lattice[LT::sup(LT::sdn(site, 0), 1)][0])*htrans(lattice[LT::sdn(site, 0)][1])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 1)][1])*(lattice[LT::sdn(LT::sdn(site, 0), 1)][0])*(lattice[LT::sdn(site, 1)][0])*(lattice[LT::sup(LT::sdn(site, 1), 0)][1])));
			long_real_t result_phi_x = real(trace(lattice[site][1]*lattice[LT::sup(site,1)][0]*htrans(lattice[LT::sup(site,0)][1])*htrans(lattice[site][0])));
			
			//long_real_t result_phi_y = real(trace((lattice[LT::sup(site, 0)][2])*htrans(lattice[LT::sup(site, 2)][0])*htrans(lattice[LT::sup(LT::sdn(site, 0), 2)][0])*htrans(lattice[LT::sdn(site, 0)][2])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 0), 2)][0])*(lattice[LT::sdn(site, 2)][0])*(lattice[LT::sup(LT::sdn(site, 2), 0)][2])));
			long_real_t result_phi_y = real(trace(lattice[site][2]*lattice[LT::sup(site,2)][0]*htrans(lattice[LT::sup(site,0)][2])*htrans(lattice[site][0])));

			//long_real_t result_phi_z = real(trace((lattice[LT::sup(site, 1)][2])*htrans(lattice[LT::sup(site, 2)][1])*htrans(lattice[LT::sup(LT::sdn(site, 1), 2)][1])*htrans(lattice[LT::sdn(site, 1)][2])*htrans(lattice[LT::sdn(LT::sdn(site, 1), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 1), 2)][1])*(lattice[LT::sdn(site, 2)][1])*(lattice[LT::sup(LT::sdn(site, 2), 1)][2])));
			long_real_t result_phi_z = real(trace(lattice[site][2]*lattice[LT::sup(site,2)][1]*htrans(lattice[LT::sup(site,1)][2])*htrans(lattice[site][1])));
			
#ifndef MULTITHREADING
			result_zero[Layout::globalIndexT(site)] += result_phi_x + result_phi_y + result_phi_z;
			result_two[Layout::globalIndexT(site)] += result_phi_x - result_phi_y;
#endif
#ifdef MULTITHREADING
			result_zero[Layout::globalIndexT(site)][omp_get_thread_num()] += result_phi_x + result_phi_y + result_phi_z;
			result_two[Layout::globalIndexT(site)][omp_get_thread_num()] += result_phi_x - result_phi_y;
#endif
		}

		//We collect the results
		for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
			zero_glueball[t] = 0;
			two_glueball[t] = 0;
			for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
				zero_glueball[t] += result_zero[t][thread];
				two_glueball[t] += result_two[t][thread];
			}
#endif
#ifndef MULTITHREADING
			zero_glueball[t] = result_zero[t];
			two_glueball[t] = result_two[t];
#endif
		}

		for (int t = 0; t < Layout::glob_t; ++t) {
			reduceAllSum(zero_glueball[t]);
			reduceAllSum(two_glueball[t]);
		}

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();

			output->push("glueball_zero");
			for (int t = 0; t < Layout::glob_t; ++t) {
				std::cout << "Glueball::O++ Operator at t " << t << " is " << zero_glueball[t]/Layout::glob_spatial_volume << std::endl;

				output->write("glueball_zero", zero_glueball[t]/Layout::glob_spatial_volume);
			}
			output->pop("glueball_zero");

			output->push("glueball_two");
			for (int t = 0; t < Layout::glob_t; ++t) {
				std::cout << "Glueball::2++ Operator at t " << t << " is " << two_glueball[t]/Layout::glob_spatial_volume << std::endl;

				output->write("glueball_two", two_glueball[t]/Layout::glob_spatial_volume);
			}
			output->pop("glueball_two");
		}
	}


	long_real_t energy_correlator[Layout::glob_t];
	long_real_t topological_correlator[Layout::glob_t];

#ifdef MULTITHREADING
	long_real_t result_energy[Layout::glob_t][omp_get_max_threads()];
	long_real_t result_topological[Layout::glob_t][omp_get_max_threads()];
	for (int i = 0; i < Layout::glob_t; ++i) {
		for (int j = 0; j < omp_get_max_threads(); ++j) {
			result_energy[i][j] = 0.;
			result_topological[i][j] = 0.;
		}
	}
#endif
#ifndef MULTITHREADING
	long_real_t result_energy[Layout::glob_t];
	long_real_t result_topological[Layout::glob_t];
	for (int i = 0; i < Layout::glob_t; ++i) {
		result_energy[i] = 0.;
		result_topological[i] = 0.;
	}
#endif

	for (int site = 0; site < lattice.localsize; ++site) {
		GaugeGroup *tmpF = new GaugeGroup[6];
		long_real_t site_energy = 0.;
		long_real_t site_topological = 0.;
		tmpF[0] = htrans(lattice[LT::sdn(site, 0)][0])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 1)][1])*(lattice[LT::sdn(LT::sdn(site, 0), 1)][0])*(lattice[LT::sdn(site, 1)][1]) + htrans(lattice[LT::sdn(site, 1)][1])*(lattice[LT::sdn(site, 1)][0])*(lattice[LT::sup(LT::sdn(site, 1), 0)][1])*htrans(lattice[site][0]) + (lattice[site][0])*(lattice[LT::sup(site, 0)][1])*htrans(lattice[LT::sup(site, 1)][0])*htrans(lattice[site][1]) + (lattice[site][1])*htrans(lattice[LT::sup(LT::sdn(site, 0), 1)][0])*htrans(lattice[LT::sdn(site, 0)][1])*(lattice[LT::sdn(site, 0)][0]);
		tmpF[1] = htrans(lattice[LT::sdn(site, 0)][0])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 0), 2)][0])*(lattice[LT::sdn(site, 2)][2]) + htrans(lattice[LT::sdn(site, 2)][2])*(lattice[LT::sdn(site, 2)][0])*(lattice[LT::sup(LT::sdn(site, 2), 0)][2])*htrans(lattice[site][0]) + (lattice[site][0])*(lattice[LT::sup(site, 0)][2])*htrans(lattice[LT::sup(site, 2)][0])*htrans(lattice[site][2]) + (lattice[site][2])*htrans(lattice[LT::sup(LT::sdn(site, 0), 2)][0])*htrans(lattice[LT::sdn(site, 0)][2])*(lattice[LT::sdn(site, 0)][0]);
		tmpF[2] = htrans(lattice[LT::sdn(site, 0)][0])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 3)][3])*(lattice[LT::sdn(LT::sdn(site, 0), 3)][0])*(lattice[LT::sdn(site, 3)][3]) + htrans(lattice[LT::sdn(site, 3)][3])*(lattice[LT::sdn(site, 3)][0])*(lattice[LT::sup(LT::sdn(site, 3), 0)][3])*htrans(lattice[site][0]) + (lattice[site][0])*(lattice[LT::sup(site, 0)][3])*htrans(lattice[LT::sup(site, 3)][0])*htrans(lattice[site][3]) + (lattice[site][3])*htrans(lattice[LT::sup(LT::sdn(site, 0), 3)][0])*htrans(lattice[LT::sdn(site, 0)][3])*(lattice[LT::sdn(site, 0)][0]);
		tmpF[3] = htrans(lattice[LT::sdn(site, 1)][1])*htrans(lattice[LT::sdn(LT::sdn(site, 1), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 1), 2)][1])*(lattice[LT::sdn(site, 2)][2]) + htrans(lattice[LT::sdn(site, 2)][2])*(lattice[LT::sdn(site, 2)][1])*(lattice[LT::sup(LT::sdn(site, 2), 1)][2])*htrans(lattice[site][1]) + (lattice[site][1])*(lattice[LT::sup(site, 1)][2])*htrans(lattice[LT::sup(site, 2)][1])*htrans(lattice[site][2]) + (lattice[site][2])*htrans(lattice[LT::sup(LT::sdn(site, 1), 2)][1])*htrans(lattice[LT::sdn(site, 1)][2])*(lattice[LT::sdn(site, 1)][1]);
		tmpF[4] = htrans(lattice[LT::sdn(site, 1)][1])*htrans(lattice[LT::sdn(LT::sdn(site, 1), 3)][3])*(lattice[LT::sdn(LT::sdn(site, 1), 3)][1])*(lattice[LT::sdn(site, 3)][3]) + htrans(lattice[LT::sdn(site, 3)][3])*(lattice[LT::sdn(site, 3)][1])*(lattice[LT::sup(LT::sdn(site, 3), 1)][3])*htrans(lattice[site][1]) + (lattice[site][1])*(lattice[LT::sup(site, 1)][3])*htrans(lattice[LT::sup(site, 3)][1])*htrans(lattice[site][3]) + (lattice[site][3])*htrans(lattice[LT::sup(LT::sdn(site, 1), 3)][1])*htrans(lattice[LT::sdn(site, 1)][3])*(lattice[LT::sdn(site, 1)][1]);
		tmpF[5] = htrans(lattice[LT::sdn(site, 2)][2])*htrans(lattice[LT::sdn(LT::sdn(site, 2), 3)][3])*(lattice[LT::sdn(LT::sdn(site, 2), 3)][2])*(lattice[LT::sdn(site, 3)][3]) + htrans(lattice[LT::sdn(site, 3)][3])*(lattice[LT::sdn(site, 3)][2])*(lattice[LT::sup(LT::sdn(site, 3), 2)][3])*htrans(lattice[site][2]) + (lattice[site][2])*(lattice[LT::sup(site, 2)][3])*htrans(lattice[LT::sup(site, 3)][2])*htrans(lattice[site][3]) + (lattice[site][3])*htrans(lattice[LT::sup(LT::sdn(site, 2), 3)][2])*htrans(lattice[LT::sdn(site, 2)][3])*(lattice[LT::sdn(site, 2)][2]);
		for (unsigned int i = 0; i < 6; ++i) {
			//Manual antialiasing, error of eigen!
			GaugeGroup antialias = tmpF[i];
			tmpF[i] = (1./8.)*(antialias - htrans(antialias));
			site_energy += real(trace(tmpF[i]*tmpF[i]));
		}
		site_topological = real(trace(tmpF[2]*tmpF[3]) - trace(tmpF[1]*tmpF[4]) + trace(tmpF[0]*tmpF[5]))/(4.*PI*PI);

#ifndef MULTITHREADING
		result_energy[Layout::globalIndexT(site)] += site_energy;
		result_topological[Layout::globalIndexT(site)] += site_topological;
#endif
#ifdef MULTITHREADING
		result_energy[Layout::globalIndexT(site)][omp_get_thread_num()] += site_energy;
		result_topological[Layout::globalIndexT(site)][omp_get_thread_num()] += site_topological;
#endif

		delete[] tmpF;
	}

	//We collect the results
	for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
		energy_correlator[t] = 0;
		topological_correlator[t] = 0;
		for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
			energy_correlator[t] += result_energy[t][thread];
			topological_correlator[t] += result_topological[t][thread];
		}
#endif
#ifndef MULTITHREADING
		energy_correlator[t] = result_energy[t];
		topological_correlator[t] = result_topological[t];
#endif
	}

	for (int t = 0; t < Layout::glob_t; ++t) {
		reduceAllSum(energy_correlator[t]);
		reduceAllSum(topological_correlator[t]);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("energy_correlator");
		output->push("topological_correlator");
		for (int t = 0; t < Layout::glob_t; ++t) {
			output->write("energy_correlator", energy_correlator[t]/Layout::glob_spatial_volume);
			output->write("topological_correlator", topological_correlator[t]/Layout::glob_spatial_volume);
		}
		output->pop("energy_correlator");
		output->pop("topological_correlator");
	}



	//then we measure the non zero projection
	bool doMomentum;
	try {
		doMomentum = (environment.configurations.get<std::string>("momentum_operator") == "true");
	} catch (NotFoundOption& err) {
		doMomentum = false;
	}

	if (doMomentum) {
		unsigned int maximum_momentum = environment.configurations.get<unsigned int>("maximum_momentum");
		for (unsigned int p_x = 0; p_x < maximum_momentum; ++p_x) {
			for (unsigned int p_y = 0; p_y < maximum_momentum; ++p_y) {
				for (unsigned int p_z = 0; p_z < maximum_momentum; ++p_z) {

					long_real_t zero_glueball[Layout::glob_t];
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
						int x = Layout::globalIndexX(site);
						int y = Layout::globalIndexY(site);
						int z = Layout::globalIndexZ(site);
						std::complex<real_t> momentum_factor(cos(2.*PI*(static_cast<real_t>(p_x*x)/Layout::glob_x + static_cast<real_t>(p_y*y)/Layout::glob_y + static_cast<real_t>(p_z*z)/Layout::glob_z)),sin(2.*PI*(static_cast<real_t>(p_x*x)/Layout::glob_x + static_cast<real_t>(p_y*y)/Layout::glob_y + static_cast<real_t>(p_z*z)/Layout::glob_z)));
						long_real_t result_phi_x = real(momentum_factor*trace((lattice[LT::sup(site, 0)][1])*htrans(lattice[LT::sup(site, 1)][0])*htrans(lattice[LT::sup(LT::sdn(site, 0), 1)][0])*htrans(lattice[LT::sdn(site, 0)][1])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 1)][1])*(lattice[LT::sdn(LT::sdn(site, 0), 1)][0])*(lattice[LT::sdn(site, 1)][0])*(lattice[LT::sup(LT::sdn(site, 1), 0)][1])));
						long_real_t result_phi_y = real(momentum_factor*trace((lattice[LT::sup(site, 0)][2])*htrans(lattice[LT::sup(site, 2)][0])*htrans(lattice[LT::sup(LT::sdn(site, 0), 2)][0])*htrans(lattice[LT::sdn(site, 0)][2])*htrans(lattice[LT::sdn(LT::sdn(site, 0), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 0), 2)][0])*(lattice[LT::sdn(site, 2)][0])*(lattice[LT::sup(LT::sdn(site, 2), 0)][2])));
						long_real_t result_phi_z = real(momentum_factor*trace((lattice[LT::sup(site, 1)][2])*htrans(lattice[LT::sup(site, 2)][1])*htrans(lattice[LT::sup(LT::sdn(site, 1), 2)][1])*htrans(lattice[LT::sdn(site, 1)][2])*htrans(lattice[LT::sdn(LT::sdn(site, 1), 2)][2])*(lattice[LT::sdn(LT::sdn(site, 1), 2)][1])*(lattice[LT::sdn(site, 2)][1])*(lattice[LT::sup(LT::sdn(site, 2), 1)][2])));
#ifndef MULTITHREADING
						result[Layout::globalIndexT(site)] += result_phi_x + result_phi_y + result_phi_z;
#endif
#ifdef MULTITHREADING
						result[Layout::globalIndexT(site)][omp_get_thread_num()] += result_phi_x + result_phi_y + result_phi_z;
#endif
					}

					//We collect the results
					for (int t = 0; t < Layout::glob_t; ++t) {
#ifdef MULTITHREADING
						zero_glueball[t] = 0;
						for (int thread = 0; thread < omp_get_max_threads(); ++thread) {
							zero_glueball[t] += result[t][thread];
						}
#endif
#ifndef MULTITHREADING
						zero_glueball[t] = result[t];
#endif
					}

					for (int t = 0; t < Layout::glob_t; ++t) {
						reduceAllSum(zero_glueball[t]);
					}

					if (environment.measurement && isOutputProcess()) {
						GlobalOutput* output = GlobalOutput::getInstance();

						output->push(std::string("glueball_zero_")+toString(p_x)+"_"+toString(p_y)+"_"+toString(p_z));
						for (int t = 0; t < Layout::pgrid_t*Layout::loc_t; ++t) {
							std::cout << "Glueball Operator at t " << t << " momentum (" << p_x <<"," << p_y << "," << p_z << ") is " << zero_glueball[t]/Layout::glob_spatial_volume << std::endl;

							output->write(std::string("glueball_zero_")+toString(p_x)+"_"+toString(p_y)+"_"+toString(p_z), zero_glueball[t]/Layout::glob_spatial_volume);
						}
						output->pop(std::string("glueball_zero_")+toString(p_x)+"_"+toString(p_y)+"_"+toString(p_z));
					}

				}
			}
		}
	}
}

} /* namespace Update */
