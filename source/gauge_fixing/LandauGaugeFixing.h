#ifndef LANDAUGAUGEFIXING_H
#define LANDAUGAUGEFIXING_H

#include "LatticeSweep.h"
#include "utils/RandomSeed.h"
#include "gauge_fixing/GaugeFixing.h"

namespace Update {

class LandauGaugeFixing : public Update::LatticeSweep, public Update::GaugeFixing {
public:
	LandauGaugeFixing();
	LandauGaugeFixing(const LandauGaugeFixing&);
	~LandauGaugeFixing();

	void execute(environment_t& environment);

	void generateOverrelaxationTransformation(extended_matrix_lattice_t& gauge_transformation, const extended_gauge_lattice_t& lattice, int d);

	static void registerParameters(std::map<std::string, Option>&);
protected:
	long_real_t functional(const extended_gauge_lattice_t& lattice);

	long_real_t deviation(const extended_gauge_lattice_t& lattice) const;

	long_real_t gaugeFixing(extended_gauge_lattice_t& lattice, const real_t& epsilon1, const real_t& beta1, const real_t& epsilon2, const real_t& beta2, const real_t& epsilon3, const real_t& beta3, unsigned int steps, unsigned int local_steps, const real_t& precision, unsigned int output_steps);
};

}

#endif

