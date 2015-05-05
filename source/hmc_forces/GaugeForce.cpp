/*
 * GaugeForce.cpp
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#include "GaugeForce.h"

namespace Update {

GaugeForce::GaugeForce() { }

GaugeForce::~GaugeForce() { }

GaugeGroup GaugeForce::force(const environment_t& environment, int site, int mu) const {
	return this->force(environment.gaugeLinkConfiguration, site, mu);
}

} /* namespace Update */
