/*
 * PolyakovLoop.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "PolyakovLoopEigenvalues.h"
#include "io/GlobalOutput.h"

namespace Update {

PolyakovLoopEigenvalues::PolyakovLoopEigenvalues() { }

PolyakovLoopEigenvalues::~PolyakovLoopEigenvalues() { }

void PolyakovLoopEigenvalues::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	typedef extended_gauge_lattice_t LT;

	extended_gauge_lattice_t tmp = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t swap;
	extended_gauge_lattice_t polyakov;

#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_identity(polyakov[site][mu]);
		}
	}

	for (int t = 0; t < Layout::glob_t; ++t) {
#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			if (Layout::globalIndexT(site) == 0) {
				polyakov[site][3] = polyakov[site][3]*tmp[site][3];
			}
		}

		//Antialias
		swap = tmp;

#pragma omp parallel for
		for (int site = 0; site < Layout::localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmp[site][mu] = swap[LT::sup(site,3)][mu];
			}
		}
		tmp.updateHalo();
	}

#ifdef MULTITHREADING
	int bins_th[315][omp_get_num_threads()];
	for (unsigned int i = 0; i < 315; ++i) {
		for (int j = 0; j < omp_get_num_threads(); ++j) {
			bins_th[i][j] = 0;
		}
	}
#endif
	int bins[315];
	for (unsigned int i = 0; i < 315; ++i) {
		bins[i] = 0;
	}

#pragma omp parallel for
	for (int site = 0; site < Layout::localsize; ++site) {
		if (Layout::globalIndexT(site) == 0) {
#ifdef EIGEN
			Eigen::ComplexEigenSolver<GaugeGroup> es(polyakov[site][3]);
			for (int i = 0; i < numberColors; ++i) {
				if (static_cast<int>(50.*(arg(es.eigenvalues()[i])+PI)) > 314 || static_cast<int>(50.*(arg(es.eigenvalues()[i])+PI)) < 0) {
					std::cout << "Fatal index error! " << static_cast<int>(50.*(arg(es.eigenvalues()[i])+PI)) << " " << arg(es.eigenvalues()[i]) << std::endl;
				}
#ifdef MULTITHREADING
				bins_th[static_cast<int>(50.*(arg(es.eigenvalues()[i])+PI))][omp_get_thread_num()] += 1;
#endif
#ifndef MULTITHREADING
				bins[static_cast<int>(50.*(arg(es.eigenvalues()[i])+PI))] += 1;
#endif
			}
#endif
#ifdef ARMADILLO
			FundamentalVector eigval;
			GaugeGroup eigvec;
			arma::eig_gen(eigval, eigvec, polyakov[site][3]);
			for (int i = 0; i < numberColors; ++i) {
				if (static_cast<int>(50.*(arg(eigval[i])+PI)) > 314 || static_cast<int>(50.*(arg(eigval[i])+PI)) < 0) {
					std::cout << "Fatal index error!" << std::endl;
				}
#ifdef MULTITHREADING
				bins_th[static_cast<int>(50.*(arg(eigval[i])+PI))][omp_get_thread_num()] += 1;
#endif
#ifndef MULTITHREADING
				bins[static_cast<int>(50.*(arg(eigval[i])+PI))] += 1;
#endif
			}
#endif
		}
	}

#ifdef MULTITHREADING
	for (unsigned int i = 0; i < 315; ++i) {
		for (int j = 0; j < omp_get_num_threads(); ++j) {
			bins[i] += bins_th[i][j];
		}
	}
#endif

	for (unsigned int i = 0; i < 315; ++i) {
		reduceAllSum(bins[i]);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("polyakov_eigenvalues");

		for (unsigned int i = 0; i < 315; ++i) {
			output->write("polyakov_eigenvalues", bins[i]);
		}

		output->pop("polyakov_eigenvalues");
	}
}

} /* namespace Update */
