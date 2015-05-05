/*
 * ReadGaugeConfiguration.cpp
 *
 *  Created on: Oct 22, 2013
 *      Author: spiem_01
 */

#include "ReadGaugeConfiguration.h"
#include "starters/ReadStartGaugeConfiguration.h"

namespace Update {

int ReadGaugeConfiguration::instanceCounter = 0;

ReadGaugeConfiguration::ReadGaugeConfiguration() : LatticeSweep(), counter(-1) {
	++instanceCounter;
	if (instanceCounter > 1 && isOutputProcess()) {
		std::cout << "ReadGaugeConfiguration::Warning, the sweep could not work meaningfully if instanced many time!" << std::endl;
	}
}

ReadGaugeConfiguration::ReadGaugeConfiguration(const ReadGaugeConfiguration& toCopy) : LatticeSweep(toCopy), counter(-1) {
	++instanceCounter;
	if (instanceCounter > 1 && isOutputProcess()) {
		std::cout << "ReadGaugeConfiguration::Warning, the sweep could not work meaningfully if instanced many time!" << std::endl;
	}
}

ReadGaugeConfiguration::~ReadGaugeConfiguration() {
	--instanceCounter;
}

void ReadGaugeConfiguration::execute(environment_t& environment) {
	if (counter == -1) counter = environment.configurations.get<unsigned int>("read_start_number");
	ReadStartGaugeConfiguration::readConfiguration(environment, counter);
	counter += environment.configurations.get<unsigned int>("read_step");
}

} /* namespace Update */
