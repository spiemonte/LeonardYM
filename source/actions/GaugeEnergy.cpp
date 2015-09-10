/*
 * GaugeEnergy.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: spiem_01
 */

#include "GaugeEnergy.h"
#include "GaugeAction.h"
#include "io/GlobalOutput.h"

namespace Update {

GaugeEnergy::GaugeEnergy() { }

GaugeEnergy::~GaugeEnergy() { }

void GaugeEnergy::execute(environment_t& environment) {
	GaugeAction* gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));
	long_real_t gaugeEnergy = gaugeAction->energy(environment);
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("gauge_energy");

		std::cout << "GaugeEnergy::Energy value " << gaugeEnergy << std::endl; //TODO
		output->write("gauge_energy", gaugeEnergy);

		output->pop("gauge_energy");
	}
	delete gaugeAction;
}

} /* namespace Update */
