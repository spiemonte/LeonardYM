#ifndef GAUGEENERGY_H_
#define GAUGEENERGY_H_

#include "LatticeSweep.h"

namespace Update {

class GaugeEnergy: public LatticeSweep {
public:
	GaugeEnergy();
	~GaugeEnergy();

	virtual void execute(environment_t& environment);
};

} /* namespace Update */
#endif /* GAUGEENERGY_H_ */
