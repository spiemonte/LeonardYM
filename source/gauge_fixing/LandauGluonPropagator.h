#ifndef LANDAUGLUONPROPAGATOR_H
#define LANDAUGLUONPROPAGATOR_H

#include "LatticeSweep.h"
#include "gauge_fixing/LandauGaugeFixing.h"

namespace Update {

class LandauGluonPropagator : public Update::LandauGaugeFixing {
public:
	LandauGluonPropagator();
	LandauGluonPropagator(const LandauGluonPropagator&);
	~LandauGluonPropagator();

	void execute(environment_t& environment);

	static void registerParameters(po::options_description&);
};

}

#endif
