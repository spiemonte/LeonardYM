#ifndef LANDAUGAUGEFIXING_H
#define LANDAUGAUGEFIXING_H

#include "LatticeSweep.h"
#include "utils/RandomSeed.h"
#include "utils/GaugeFixing.h"

namespace Update {

class LandauGaugeFixing : public Update::LatticeSweep, public Update::GaugeFixing {
public:
	LandauGaugeFixing();
	LandauGaugeFixing(const LandauGaugeFixing&);
	~LandauGaugeFixing();

	void execute(environment_t& environment);

	void generateOverrelaxationTransformation(extended_matrix_lattice_t& gauge_transformation, const extended_gauge_lattice_t& lattice, int d);

	static void registerParameters(po::options_description&);
protected:
	long_real_t functional(const extended_gauge_lattice_t& lattice);
};

}

#endif

