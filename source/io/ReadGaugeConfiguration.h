/*
 * ReadGaugeConfiguration.h
 *
 *  Created on: Oct 22, 2013
 *      Author: spiem_01
 */

#ifndef READGAUGECONFIGURATION_H_
#define READGAUGECONFIGURATION_H_
#include "LatticeSweep.h"

namespace Update {

class ReadGaugeConfiguration : public LatticeSweep {
public:
	ReadGaugeConfiguration();
	ReadGaugeConfiguration(const ReadGaugeConfiguration&);
	virtual ~ReadGaugeConfiguration();

	virtual void execute(environment_t& environment);

private:
	static int instanceCounter;
	int counter;
};

} /* namespace Update */
#endif /* READGAUGECONFIGURATION_H_ */
