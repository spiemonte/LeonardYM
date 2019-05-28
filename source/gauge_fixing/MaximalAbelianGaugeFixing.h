#ifndef MAXIMALABELIANGAUGEFIXING_H
#define MAXIMALABELIANGAUGEFIXING_H

#include "LatticeSweep.h"
#include "utils/RandomSeed.h"
#include "gauge_fixing/GaugeFixing.h"

namespace Update {

class MaximalAbelianGaugeFixing : public Update::LatticeSweep, public Update::GaugeFixing {
public:
	MaximalAbelianGaugeFixing();
	MaximalAbelianGaugeFixing(const MaximalAbelianGaugeFixing&);
	~MaximalAbelianGaugeFixing();

	void execute(environment_t& environment);

	static void registerParameters(po::options_description&);
protected:
	long_real_t functional(const extended_gauge_lattice_t& lattice);
};

}

#endif

