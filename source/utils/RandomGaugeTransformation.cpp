/*
 * RandomStartGaugeConfiguration.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#include "RandomGaugeTransformation.h"
#include "MatrixTypedef.h"
#include "utils/ReUnit.h"

namespace Update {

RandomGaugeTransformation::RandomGaugeTransformation() : LatticeSweep(), rng(RandomSeed::randomSeed()),
#if NUMCOLORS > 2
		randomNormal(RandomSeed::getNormalNumberGenerator(rng))
#endif
#if NUMCOLORS == 2
		randomUniform(RandomSeed::getRandomNumberGenerator(rng))
#endif
	{ }

RandomGaugeTransformation::~RandomGaugeTransformation() { }

void RandomGaugeTransformation::execute(environment_t& environment) {
	//The random gauge transformation
	extended_matrix_lattice_t omega;
	typedef extended_matrix_lattice_t LT;

	//First generate the gauge transformation
	for (int site = 0; site < omega.localsize; ++site) {
#if NUMCOLORS > 2
#ifdef EIGEN
		GaugeGroup random;
		//Fill the matrix with element from gaussian distribution (mu = 0, sigma = 1)
		for (int i = 0; i < numberColors; ++i) {
			for (int j = 0; j < numberColors; ++j) {
				random(i,j) = std::complex<real_t>(randomNormal(),randomNormal()); //call the operator()
			}
		}
		Eigen::HouseholderQR<GaugeGroup> qrdecomp(random);
		omega[site] = qrdecomp.householderQ();
		//Now set the det(Q) to 1:
		complex determinant = conj(det(omega[site]));
		// by multiply the first row with the complex conjugate of det:
		for (int i = 0; i < numberColors; ++i) {
			omega[site].at(0,i) *= determinant;
		}
#endif
#ifdef ARMADILLO
		//Fill the matrix with element from gaussian distribution (mu = 0, sigma = 1)
		for (int i = 0; i < numberColors; ++i) {
			for (int j = 0; j < numberColors; ++j) {
				omega(i,j) = std::complex<real_t>(randomNormal(),randomNormal()); //call the operator()
			}
		}
		//TODO: not really working for SU(N)
#endif
#endif
#if NUMCOLORS == 2
		real_t chi12 = (2*PI*randomUniform());
		real_t psi12 = (2*PI*randomUniform());
		real_t phi12 = asin(pow(randomUniform(), (1./2.)));
		omega[site].at(0,0) = (static_cast<real_t>(cos(phi12))*complex(cos(psi12), sin(psi12)));
		omega[site].at(0,1) = (complex(cos(chi12), sin(chi12))*static_cast<real_t>(sin(phi12)));
		omega[site].at(1,0) = (-(complex(cos(chi12), (-sin(chi12)))*static_cast<real_t>(sin(phi12))));
		omega[site].at(1,1) = (static_cast<real_t>(cos(phi12))*complex(cos(psi12), (-sin(psi12))));
#endif
	}
	omega.updateHalo();

	//Perform the gauge transformation
#pragma omp parallel for
	for (int site = 0; site < omega.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup swap = environment.gaugeLinkConfiguration[site][mu];
			environment.gaugeLinkConfiguration[site][mu] = omega[site]*swap*htrans(omega[LT::sup(site,mu)]);
		}
	}

	environment.gaugeLinkConfiguration.updateHalo();

	//Synchronize the eventual adjoint lattice
	environment.synchronize();
}

} /* namespace Update */
