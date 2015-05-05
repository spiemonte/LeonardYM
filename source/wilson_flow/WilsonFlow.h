/*
 * WilsonFlow.h
 *
 *  Created on: Oct 23, 2013
 *      Author: spiem_01
 */

#ifndef WILSONFLOW_H_
#define WILSONFLOW_H_
#include "LatticeSweep.h"
#include "actions/GaugeAction.h"
#include <utility>

namespace Update {

class WilsonFlow : public LatticeSweep {
public:
	WilsonFlow();
	WilsonFlow(const WilsonFlow&);
	virtual ~WilsonFlow();


	virtual void execute(environment_t& environment);
protected:
	void integrate(const extended_gauge_lattice_t& initialLattice, extended_gauge_lattice_t& finalLattice, GaugeAction* action, real_t time, int nSteps);
	void getForce(const extended_gauge_lattice_t& lattice, extended_gauge_lattice_t& force, GaugeAction* action);
	GaugeGroup exponential(const GaugeGroup& link, const GaugeGroup& force, real_t epsilon);

	void measureEnergy(const extended_gauge_lattice_t& lattice);

	void threeDimensionalEnergyTopologicalPlot(const extended_gauge_lattice_t& lattice, environment_t& environment);

	std::pair<long_real_t,long_real_t> measureEnergyAndTopologicalCharge(const extended_gauge_lattice_t& lattice, int site);

private:
	real_t gaugeEnergy;
	real_t topologicalCharge;

	long_real_t *energy_correlator;
	long_real_t *topological_correlator;
};

} /* namespace Update */
#endif /* WILSONFLOW_H_ */
