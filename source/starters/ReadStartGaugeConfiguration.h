/*
 * ReadStartGaugeConfiguration.h
 *
 *  Created on: May 28, 2012
 *      Author: spiem_01
 */

#ifndef READSTARTGAUGECONFIGURATION_H_
#define READSTARTGAUGECONFIGURATION_H_

#include "StartGaugeConfiguration.h"

namespace Update {

class ReadStartGaugeConfiguration : public Update::StartGaugeConfiguration {
public:
	ReadStartGaugeConfiguration();
	~ReadStartGaugeConfiguration();

	/**
	 * This function initialize the fundamental gauge link configuration accordingly to that stored in file "start_gauge_configuration_file"
	 * @param enviroment
	 * @param sweep
	 * @param n
	 */
	virtual void execute(environment_t& environment);

	static bool readConfiguration(environment_t& environment, int numberConfiguration);

private:
	/*
#ifndef TESTLATTICE
	static bool isreadingrank_() {
		return (Math::LinkConf::LinkConfiguration<GaugeGroup>::Layout::isIoRank());
	}
#endif
	*/
};

} /* namespace Update */
#endif /* READSTARTGAUGECONFIGURATION_H_ */
