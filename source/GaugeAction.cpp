/*
 * GaugeAction.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#include "GaugeAction.h"
#include "WilsonGaugeAction.h"
#include "ImprovedGaugeAction.h"

namespace Update {

GaugeAction::GaugeAction(double _beta) : Energy(), beta(_beta) { }

GaugeAction::~GaugeAction() { }

GaugeAction* GaugeAction::getInstance(const std::string& name, double _beta) {
	if (name == "StandardWilson") {
		return new WilsonGaugeAction(_beta);
	} else if (name == "Improved") {
		return new ImprovedGaugeAction(_beta);
	} else {
		if (isOutputProcess()) std::cout << "Action name unknown: " << name << std::endl;
		exit(1);
	}
}

void GaugeAction::setBeta(real_t _beta) {
	beta = _beta;
}

real_t GaugeAction::getBeta() const {
	return beta;
}

} /* namespace Update */
