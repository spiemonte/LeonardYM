/*
 * GaugeAction.h
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifndef GAUGEACTION_H_
#define GAUGEACTION_H_

#include "Environment.h"
#include "GaugeForce.h"
#include "Energy.h"
#include <string>

namespace Update {

class GaugeAction : public Energy, public GaugeForce {
public:
	GaugeAction(double _beta);
	virtual ~GaugeAction();

	static GaugeAction* getInstance(const std::string& name, double beta);

	virtual GaugeGroup staple(const extended_gauge_lattice_t& lattice, int site, int mu) const = 0;

	virtual real_t deltaAction(const extended_gauge_lattice_t& lattice, const GaugeGroup& trial, const GaugeGroup& staple, int site, int mu) const = 0;

	void setBeta(real_t _beta);
	real_t getBeta() const;
private:
	real_t beta;
};

} /* namespace Update */
#endif /* GAUGEACTION_H_ */
