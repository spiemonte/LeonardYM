/*
 * RandomStartGaugeConfiguration.h
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifndef COLDSTARTGAUGECONFIGURATION_H_
#define COLDSTARTGAUGECONFIGURATION_H_

#include "StartGaugeConfiguration.h"

namespace Update {

class ColdStartGaugeConfiguration : public Update::StartGaugeConfiguration {
public:
	ColdStartGaugeConfiguration();
	~ColdStartGaugeConfiguration();

	/**
	 * This function initialize the fundamental gauge link configuration with a cold start, every matrix is set to unity
	 * @param enviroment
	 * @param sweep
	 * @param n
	 */
	virtual void execute(environment_t& environment);

};

} /* namespace Update */
#endif /* COLDSTARTGAUGECONFIGURATION_H_ */
