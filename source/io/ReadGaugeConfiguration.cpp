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
	bool success = ReadStartGaugeConfiguration::readConfiguration(environment, counter);
	int reattempt = 1;
	while (!success) {
		if (isOutputProcess()) std::cout << "ReadGaugeConfiguration::Warning, configuration number " << counter << " not readable!" << std::endl;
		counter += environment.configurations.get<unsigned int>("read_step");
		success = ReadStartGaugeConfiguration::readConfiguration(environment, counter);
		++reattempt;
		if (reattempt == 100) {
			if (isOutputProcess()) std::cout << "ReadGaugeConfiguration::Reached maximum reattempt iterations!" << std::endl;
			exit(73);
			break;
		}
	}
	counter += environment.configurations.get<unsigned int>("read_step");
}

} /* namespace Update */
