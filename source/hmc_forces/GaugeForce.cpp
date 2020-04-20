#include "GaugeForce.h"

namespace Update {

GaugeForce::GaugeForce() { }

GaugeForce::~GaugeForce() { }

GaugeGroup GaugeForce::force(const environment_t& environment, int site, int mu) const {
	return this->force(environment.gaugeLinkConfiguration, site, mu);
}

} /* namespace Update */
