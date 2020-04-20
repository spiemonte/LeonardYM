#ifndef STARTGAUGECONFIGURATION_H_
#define STARTGAUGECONFIGURATION_H_

#include "LatticeSweep.h"

namespace Update {

class StartGaugeConfiguration: public Update::LatticeSweep {
public:
	StartGaugeConfiguration();
	~StartGaugeConfiguration();

	/**
	 * This factory returns back an instance to the sweep that initialize the fundamental lattice
	 * @param name
	 * @return
	 */
	static StartGaugeConfiguration* getInstance(const std::string& name);
};

} /* namespace Update */
#endif /* STARTGAUGECONFIGURATION_H_ */
