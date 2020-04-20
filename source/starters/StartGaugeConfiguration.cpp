#include "StartGaugeConfiguration.h"
#include "HotStartGaugeConfiguration.h"
#include "ColdStartGaugeConfiguration.h"
#include "ReadStartGaugeConfiguration.h"

namespace Update {

StartGaugeConfiguration::StartGaugeConfiguration() { }

StartGaugeConfiguration::~StartGaugeConfiguration() { }

StartGaugeConfiguration* StartGaugeConfiguration::getInstance(const std::string& name) {
	if (name == "hotstart") {
		return new HotStartGaugeConfiguration();
	} else if (name == "coldstart") {
		return new ColdStartGaugeConfiguration();
	} else if (name == "readstart") {
		return new ReadStartGaugeConfiguration();
	} else {
		std::cout << "Unknown starter method name: " << name << std::endl;
		exit(1);
	}
}

} /* namespace Update */
