/*
 * HotStartGaugeConfiguration.h
 *
 *  Created on: Feb 28, 2012
 *      Author: spiem_01
 */

#ifndef HOTSTARTGAUGECONFIGURATION_H_
#define HOTSTARTGAUGECONFIGURATION_H_

#include "StartGaugeConfiguration.h"
#include "utils/RandomSeed.h"

namespace Update {

class HotStartGaugeConfiguration: public Update::StartGaugeConfiguration {
public:
	HotStartGaugeConfiguration();
	~HotStartGaugeConfiguration();

	/**
	 * This function initialize the fundamental gauge link configuration with an hot-start
	 * @param enviroment
	 * @param sweep
	 * @param n
	 */
	virtual void execute(environment_t& environment);
private:
	//The random number generator
	random_generator_t rng;
#if NUMCOLORS > 2
	//The random normal number generator, needed by SUN
	random_normal_generator_t randomNormal;
#endif
#if NUMCOLORS == 2
	//The random uniform generator in [0,1], needed by SU2
	random_uniform_generator_t randomUniform;
#endif
};

} /* namespace Update */
#endif /* HOTSTARTGAUGECONFIGURATION_H_ */
