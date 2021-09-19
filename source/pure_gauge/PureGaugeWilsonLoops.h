#ifndef PUREGAUGEWILSONLOOPS_H_
#define PUREGAUGEWILSONLOOPS_H_

#include "LatticeSweep.h"
#include "actions/GaugeAction.h"
#include "PureGaugeUpdater.h"
#include "PureGaugeOverrelaxation.h"

namespace Update {

class SpatialWilsonLoop {
public:
	SpatialWilsonLoop(int _x0, int _y0, int _z0, int _t0, int _R, int _T) :
		x0(_x0), y0(_y0), z0(_z0), t0(_t0), R(_R), T(_T) { }

	GaugeGroup measure(const extended_gauge_lattice_t& lattice) {
		GaugeGroup result;
		set_to_identity(result);
		typedef extended_gauge_lattice_t::Layout Layout;
		int site = Layout::getGlobalCoordinate(x0,y0,z0,t0);
#ifdef ENABLE_MPI
		exit(255); //Working only in multithreading mode
#endif
		for (int dx = 0; dx < R; ++dx) {
			result *= lattice[site][0];
			site = extended_gauge_lattice_t::sup(site,0);
		}
		for (int dy = 0; dy < T; ++dy) {
			result *= lattice[site][1];
			site = extended_gauge_lattice_t::sup(site,1);
		}
		for (int dx = R - 1; dx >= 0; --dx) {
			site = extended_gauge_lattice_t::sdn(site,0);
			result *= htrans(lattice[site][0]);
		}
		for (int dy = T - 1; dy >= 0; --dy) {
			site = extended_gauge_lattice_t::sdn(site,1);
			result *= htrans(lattice[site][1]);
		}
		return result;
	}

	int getR() const { return R; }
	int getT() const { return T; }

private:
	int x0;
	int y0;
	int z0;
	int t0;
	int R;
	int T;
};

class PureGaugeWilsonLoops : public Update::LatticeSweep {
public:
	PureGaugeWilsonLoops();
	PureGaugeWilsonLoops(const PureGaugeWilsonLoops& toCopy);
	~PureGaugeWilsonLoops();

	virtual void execute(environment_t& environment);

	static void registerParameters(std::map<std::string, Option>& desc);

protected:
	void updateSlices(environment_t& environment, GaugeAction* action, int sliceSize);

private:
	std::vector<SpatialWilsonLoop>* wilsonLoops;
	real_t* results;

	PureGaugeUpdater* pureGaugeUpdater;
	PureGaugeOverrelaxation* pureGaugeOverrelaxation;
};

} /* namespace Update */
#endif /* PUREGAUGEWILSONLOOPS_H_ */
