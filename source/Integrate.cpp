/*
 * Integrate.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#include "Integrate.h"
#include "LeapFrog.h"
#include "FourthOrderLeapFrog.h"
#include "SixthOrderLeapFrog.h"
#include "OmelyanLeapFrog.h"
#include "FourthOmelyanLeapFrog.h"
#ifdef EIGEN
#include <Eigen/Eigenvalues>
#endif
#ifdef ARMADILLO

#endif
#include "ToString.h"

namespace Update {

Integrate::Integrate() { }

Integrate::~Integrate() { }

void Integrate::updateLinkConfiguration(extended_gauge_lattice_t& linkConfiguration, const extended_gauge_lattice_t& momenta, real_t epsilon) {
#pragma omp parallel for
	for (int site = 0; site < linkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef EIGEN
			Eigen::ComplexEigenSolver<GaugeGroup> es(momenta[site][mu]);
			GaugeGroup update;
			update.zeros();
			for (int i = 0; i < numberColors; ++i) {
				update.at(i,i) = exp(epsilon*es.eigenvalues()[i]);
			}
			GaugeGroup m = es.eigenvectors();
			GaugeGroup updatenew = m * update * htrans(m);
#endif
#ifdef ARMADILLO
			FundamentalVector eigval;
			GaugeGroup eigvec;
			arma::eig_gen(eigval, eigvec, momenta[site][mu]);
			GaugeGroup update;
			set_to_zero(update);
			for (unsigned int i = 0; i < numberColors; ++i) {
				update.at(i,i) = exp(epsilon*eigval[i]);
			}
			GaugeGroup updatenew = eigvec * update * htrans(eigvec);
#endif
#ifdef MATRIXTOOLKIT
			FundamentalVector eigval;
			GaugeGroup eigvec;
			matrix_toolkit::eigensystem(eigval, eigvec, momenta[site][mu]);
			GaugeGroup update;
			set_to_zero(update);
			for (unsigned int i = 0; i < numberColors; ++i) {
				update.at(i,i) = exp(epsilon*eigval[i]);
			}
			GaugeGroup updatenew = eigvec * update * htrans(eigvec);
#endif
			linkConfiguration[site][mu] = updatenew*linkConfiguration[site][mu];
		}
	}
	//linkConfiguration.updateHalo(); Synchronize to be called outside
}

Integrate* Integrate::getInstance(const std::string & nameAlgorithm) {
	if (nameAlgorithm == "second_order") {
		return new LeapFrog();
	} else if (nameAlgorithm == "fourth_order") {
		return new FourthOrderLeapFrog();
	} else if (nameAlgorithm == "sixth_order") {
		return new SixthOrderLeapFrog();
	} else if (nameAlgorithm == "omelyan") {
		return new OmelyanLeapFrog();
	} else if (nameAlgorithm == "fourth_omelyan") {
		return new FourthOmelyanLeapFrog();
	} else {
		if (isOutputProcess()) std::cout << "Unknown integrate algorithm name " << nameAlgorithm << std::endl;
		exit(1);
	}
}

void Integrate::updateMomenta(extended_gauge_lattice_t& momenta, const extended_gauge_lattice_t& force, real_t epsilon) {
#pragma omp parallel for
	for (int site = 0; site < momenta.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			momenta[site][mu] = momenta[site][mu] - epsilon*force[site][mu];
		}
	}
	//momenta.updateHalo(); not needed anymore
}

} /* namespace Update */
