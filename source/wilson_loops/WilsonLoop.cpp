#include "WilsonLoop.h"
#include "io/GlobalOutput.h"
#include "utils/StoutSmearing.h"

namespace Update {

WilsonLoop::WilsonLoop() : tWilsonLine(0), xWilsonLine(0) { }

WilsonLoop::WilsonLoop(const WilsonLoop& toCopy) : LatticeSweep(toCopy), tWilsonLine(0), xWilsonLine(0) { }

WilsonLoop::~WilsonLoop() {
	if (tWilsonLine != 0) {
		for (int i = 0; i < TMax; ++i) {
			delete[] tWilsonLine[i];
		}
		delete[] tWilsonLine;
	}
	if (xWilsonLine != 0) {
		for (int i = 0; i < RMax; ++i) {
			delete[] xWilsonLine[i];
		}
		delete[] xWilsonLine;
	}
}

void WilsonLoop::execute(environment_t& environment) {
	//We work with reduced halos
	reduced_gauge_lattice_t originalLattice = environment.gaugeLinkConfiguration;

	unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("WilsonLoop::level_stout_smearing");
	if (numberLevelSmearing > 0) {
		double smearingRho = environment.configurations.get<double>("WilsonLoop::rho_stout_smearing");
		extended_gauge_lattice_t smearedConfiguration;
		StoutSmearing stoutSmearing;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, smearedConfiguration, numberLevelSmearing, smearingRho);
		originalLattice = smearedConfiguration;
	} else {
		if (isOutputProcess()) std::cout << "WilsonLoop::No smearing!" << std::endl;
	}

	RMax = environment.configurations.get<unsigned int>("WilsonLoop::max_r");
	TMax = environment.configurations.get<unsigned int>("WilsonLoop::max_t");

	int t_dir = environment.configurations.get<unsigned int>("WilsonLoop::t_dir");

	//tWilsonLine[T][R] contains the wilson line in the t-direction long T and shifted of R sites in the x-direction
	if (tWilsonLine == 0) {
		tWilsonLine = new reduced_matrix_lattice_t*[TMax];
		for (int i = 0; i < TMax; ++i) {
			tWilsonLine[i] = new reduced_matrix_lattice_t[RMax+1];
		}
	}
	//xWilsonLine[R][T] contains the wilson line in the x-direction long R and shifted of T sites in the t-direction
	if (xWilsonLine == 0){
		xWilsonLine = new reduced_matrix_lattice_t*[RMax];
		for (int i = 0; i < RMax; ++i) {
			xWilsonLine[i] = new reduced_matrix_lattice_t[TMax+1];
		}
	}

	typedef reduced_matrix_lattice_t LT;

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("wilson_loops");
	}

#pragma omp parallel for
	for (int site = 0; site < originalLattice.localsize; ++site) {
		tWilsonLine[0][0][site] =  originalLattice[site][t_dir];
		xWilsonLine[0][0][site] =  originalLattice[site][0];
	}
	tWilsonLine[0][0].updateHalo();
	xWilsonLine[0][0].updateHalo();

	reduced_matrix_lattice_t shiftedLattice = tWilsonLine[0][0], swap;
	for (int t = 1; t < TMax; ++t) {
#pragma omp parallel for
		for (int site = 0; site < originalLattice.localsize; ++site) {
			tWilsonLine[t][0][site] = tWilsonLine[t-1][0][site]*shiftedLattice[LT::sup(site,t_dir)];
			swap[site] = shiftedLattice[LT::sup(site,t_dir)];
		}
		swap.updateHalo();
		shiftedLattice = swap;
	}

	shiftedLattice = xWilsonLine[0][0];
	for (int r = 1; r < RMax; ++r) {
#pragma omp parallel for
		for (int site = 0; site < originalLattice.localsize; ++site) {
			xWilsonLine[r][0][site] = xWilsonLine[r-1][0][site]*shiftedLattice[LT::sup(site,0)];
			swap[site] = shiftedLattice[LT::sup(site,0)];
		}
		swap.updateHalo();
		shiftedLattice = swap;
	}

	for (int r = 0; r < RMax; ++r) {
		for (int t = 1; t < TMax + 1; ++t) {
#pragma omp parallel for
			for (int site = 0; site < originalLattice.localsize; ++site) {
				xWilsonLine[r][t][site] = xWilsonLine[r][t-1][LT::sup(site,t_dir)];
			}
			xWilsonLine[r][t].updateHalo();
		}
	}

	for (int t = 0; t < TMax; ++t) {
		for (int r = 1; r < RMax + 1; ++r) {
#pragma omp parallel for
			for (int site = 0; site < originalLattice.localsize; ++site) {
				tWilsonLine[t][r][site] = tWilsonLine[t][r-1][LT::sup(site,0)];
			}
			tWilsonLine[t][r].updateHalo();
		}
	}

	//Now we can compute the wilson loops
	for (int R = 0; R < RMax; ++R) {
		for (int T = 0; T < TMax; ++T) {
			long_real_t result = 0.;
#pragma omp parallel for reduction(+:result)
			for (int site = 0; site < originalLattice.localsize; ++site) {
				result += real(trace(xWilsonLine[R][0][site]*tWilsonLine[T][R+1][site]*htrans(xWilsonLine[R][T+1][site])*htrans(tWilsonLine[T][0][site])));
			}
			reduceAllSum(result);

			result = result/static_cast<long_real_t>(numberColors*reduced_matrix_lattice_t::Layout::globalVolume);

			if (environment.measurement && isOutputProcess()) {
				GlobalOutput* output = GlobalOutput::getInstance();

				std::cout << "Temporal Wilson loop (" << R + 1 << "," << T + 1<< "): " << result << std::endl;
				output->push("wilson_loops");
				output->write("wilson_loops", R + 1);
				output->write("wilson_loops", T + 1);
				output->write("wilson_loops", result);
				output->pop("wilson_loops");
			}
		}
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->pop("wilson_loops");
	}
}

void WilsonLoop::registerParameters(std::map<std::string, Option>& desc) {
	desc["WilsonLoop::max_t"] = Option("WilsonLoop::max_t", 6, "The maximum dimension of the Wilson loop in the R direction");
	desc["WilsonLoop::max_r"] = Option("WilsonLoop::max_r", 6, "The maximum dimension of the Wilson loop in the R direction");
	desc["WilsonLoop::t_dir"] = Option("WilsonLoop::t_dir", 3, "The time direction of the loops");
	desc["WilsonLoop::level_stout_smearing"] = Option("WilsonLoop::level_stout_smearing", 10, "Number of levels of the stout smearing");
	desc["WilsonLoop::rho_stout_smearing"] = Option("WilsonLoop::rho_stout_smearing", 0.05, "Rho stout smearing");
}

} /* namespace Update */
