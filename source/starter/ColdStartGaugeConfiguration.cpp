/*
 * RandomStartGaugeConfiguration.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#include "ColdStartGaugeConfiguration.h"
#include "MatrixTypedef.h"

namespace Update {

ColdStartGaugeConfiguration::ColdStartGaugeConfiguration() : StartGaugeConfiguration() { }

ColdStartGaugeConfiguration::~ColdStartGaugeConfiguration() { }

void ColdStartGaugeConfiguration::execute(environment_t& environment) {
#pragma omp parallel for
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			//Set to unit matrix
			for (int i = 0; i < numberColors; ++i) {
				for (int j = 0; j < numberColors; ++j) {
					environment.gaugeLinkConfiguration[site][mu].at(i,j) = std::complex<real_t>(0.,0.);
				}
				environment.gaugeLinkConfiguration[site][mu].at(i,i) = std::complex<real_t>(1.,0.);
			}
		}
	}
	environment.gaugeLinkConfiguration.updateHalo();
	//Synchronize the eventual adjoint lattice
	environment.synchronize();
}

} /* namespace Update */
