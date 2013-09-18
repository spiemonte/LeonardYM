/*
 * GaugeForce.h
 *
 *  Created on: Mar 21, 2012
 *      Author: spiem_01
 */

#ifndef GAUGEFORCE_H_
#define GAUGEFORCE_H_

#include "Force.h"

namespace Update {

class GaugeForce : public Update::Force {
public:
	GaugeForce();
	virtual ~GaugeForce();

	virtual GaugeGroup force(const environment_t& environment, int site, int mu) const;

	virtual GaugeGroup force(const extended_gauge_lattice_t& lattice, int site, int mu) const = 0;
};

} /* namespace Update */
#endif /* GAUGEFORCE_H_ */
