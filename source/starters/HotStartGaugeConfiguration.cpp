/*
 * RandomStartGaugeConfiguration.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#include "HotStartGaugeConfiguration.h"
#include "MatrixTypedef.h"
#include "utils/ReUnit.h"
#ifndef PI
#define PI 3.14159265358979323846264338327950288419
#endif

namespace Update {

HotStartGaugeConfiguration::HotStartGaugeConfiguration() : StartGaugeConfiguration(), rng(RandomSeed::randomSeed()),
#if NUMCOLORS > 2
		randomNormal(RandomSeed::getNormalNumberGenerator(rng))
#endif
#if NUMCOLORS == 2
		randomUniform(RandomSeed::getRandomNumberGenerator(rng))
#endif
	{ }

HotStartGaugeConfiguration::~HotStartGaugeConfiguration() { }

void HotStartGaugeConfiguration::execute(environment_t& environment) {
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
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
			environment.gaugeLinkConfiguration[site][mu] = qrdecomp.householderQ();
			//Now set the det(Q) to 1:
			complex determinant = conj(det(environment.gaugeLinkConfiguration[site][mu]));
			// by multiply the first row with the complex conjugate of det:
			for (int i = 0; i < numberColors; ++i) {
				environment.gaugeLinkConfiguration[site][mu].at(0,i) *= determinant;
			}
#endif
#ifdef ARMADILLO
			GaugeGroup random;
			//Fill the matrix with element from gaussian distribution (mu = 0, sigma = 1)
			for (int i = 0; i < numberColors; ++i) {
				for (int j = 0; j < numberColors; ++j) {
					random(i,j) = std::complex<real_t>(randomNormal(),randomNormal()); //call the operator()
				}
			}
#endif
#endif
#if NUMCOLORS == 2
			real_t chi12 = (2*PI*randomUniform());
			real_t psi12 = (2*PI*randomUniform());
			real_t phi12 = asin(pow(randomUniform(), (1./2.)));
			environment.gaugeLinkConfiguration[site][mu].at(0,0) = (static_cast<real_t>(cos(phi12))*complex(cos(psi12), sin(psi12)));
			environment.gaugeLinkConfiguration[site][mu].at(0,1) = (complex(cos(chi12), sin(chi12))*static_cast<real_t>(sin(phi12)));
			environment.gaugeLinkConfiguration[site][mu].at(1,0) = (-(complex(cos(chi12), (-sin(chi12)))*static_cast<real_t>(sin(phi12))));
			environment.gaugeLinkConfiguration[site][mu].at(1,1) = (static_cast<real_t>(cos(phi12))*complex(cos(psi12), (-sin(psi12))));
#endif
		}
	}
	environment.gaugeLinkConfiguration.updateHalo();
	//Synchronize the eventual adjoint lattice
	environment.synchronize();
	ReUnit reunit;
	reunit.execute(environment);
}

} /* namespace Update */
