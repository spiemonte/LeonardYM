/*
 * StoutSmearing.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#ifndef STOUTSMEARING_H_
#define STOUTSMEARING_H_
#include "Environment.h"

namespace Update {

class StoutSmearing {
public:
	StoutSmearing();
	~StoutSmearing();

	void smearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, real_t rho);
	void spatialSmearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, unsigned int numberLevels, real_t rho);
private:
	GaugeGroup exp(const GaugeGroup& toExp);
};

} /* namespace Update */
#endif /* STOUTSMEARING_H_ */
