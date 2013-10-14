/*
 * PureGaugeWilsonLoops.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: spiem_01
 */

#include "PureGaugeWilsonLoops.h"
#include "GlobalOutput.h"
#include "Checkerboard.h"

namespace Update {

class TwoLinkOperator {
	public:
		TwoLinkOperator() { }

		void setToZero() {
			for (int a = 0; a < numberColors; ++a) {
				for (int b = 0; b < numberColors; ++b) {
					for (int c = 0; c < numberColors; ++c) {
						for (int d = 0; d < numberColors; ++d) {
							twoLink[a][b][c][d] = 0.;
						}
					}
				}
			}
		}

		void add(const GaugeGroup& U1, const GaugeGroup& U2) {
			GaugeGroup hU2 = htrans(U2);
			for (int a = 0; a < numberColors; ++a) {
				for (int b = 0; b < numberColors; ++b) {
					for (int c = 0; c < numberColors; ++c) {
						for (int d = 0; d < numberColors; ++d) {
							twoLink[a][b][c][d] += U1.at(a,b)*hU2.at(c,d);
						}
					}
				}
			}
		}

		const std::complex<real_t>& at(int a, int b, int c, int d) const {
			return twoLink[a][b][c][d];
		}

		std::complex<real_t>& at(int a, int b, int c, int d) {
			return twoLink[a][b][c][d];
		}

		TwoLinkOperator operator*(const TwoLinkOperator& snd) {
			TwoLinkOperator result;
			for (int a = 0; a < numberColors; ++a) {
				for (int b = 0; b < numberColors; ++b) {
					for (int c = 0; c < numberColors; ++c) {
						for (int d = 0; d < numberColors; ++d) {
							result.at(a,b,c,d) = 0.;
							for (int e = 0; e < numberColors; ++e) {
								for (int f = 0; f < numberColors; ++f) {
									result.at(a,b,c,d) += this->at(a,e,c,f)*snd.at(e,b,f,d);
								}
							}
						}
					}
				}
			}
			return result;	
		}

		void normalize(int factor) {
			for (int a = 0; a < numberColors; ++a) {
				for (int b = 0; b < numberColors; ++b) {
					for (int c = 0; c < numberColors; ++c) {
						for (int d = 0; d < numberColors; ++d) {
							twoLink[a][b][c][d] /= factor;
						}
					}
				}
			}
		}
	private:
		std::complex<real_t> twoLink[numberColors][numberColors][numberColors][numberColors];
};

PureGaugeWilsonLoops::PureGaugeWilsonLoops() : LatticeSweep(), wilsonLoops(0), results(0), pureGaugeUpdater(new PureGaugeUpdater()), pureGaugeOverrelaxation(new PureGaugeOverrelaxation()) { }

PureGaugeWilsonLoops::PureGaugeWilsonLoops(const PureGaugeWilsonLoops& toCopy) : LatticeSweep(toCopy), wilsonLoops(0), results(0), pureGaugeUpdater(new PureGaugeUpdater()), pureGaugeOverrelaxation(new PureGaugeOverrelaxation()) { }

PureGaugeWilsonLoops::~PureGaugeWilsonLoops() {
	if (results) delete[] results;
	if (wilsonLoops) delete[] wilsonLoops;
	delete pureGaugeUpdater;
	delete pureGaugeOverrelaxation;
}

void PureGaugeWilsonLoops::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout LT;

#ifndef ENABLE_MPI

	int RMax = environment.configurations.get<unsigned int>("max_r_wilsonloop");
	int RMin = environment.configurations.get<unsigned int>("min_r_wilsonloop");
	int TMax = environment.configurations.get<unsigned int>("max_t_wilsonloop");
	int TMin = environment.configurations.get<unsigned int>("min_t_wilsonloop");

	//Collect together all the wilson loop that should be measured
	if (wilsonLoops == 0) {
		wilsonLoops = new std::vector<SpatialWilsonLoop>[(RMax-RMin+1)*(TMax-TMin+1)];
		results = new real_t[(RMax-RMin+1)*(TMax-TMin+1)];
		int index = 0;
		for (int R = RMin; R <= RMax; ++R) {
			for (int T = TMin; T <= TMax; ++T) {
				for (int x = 0; x < LT::glob_x; x += 2) {
					for (int y = 0; y < LT::glob_y; y += 2) {
						for (int z = 0; z < LT::glob_z; z += 2) {
							for (int t = 0; t < LT::glob_t; t += 2) {
								wilsonLoops[index].push_back(SpatialWilsonLoop(x,y,z,t,R,T));
							}
						}
					}
				}
				++index;
			}
		}
	}

	//measure them
	for (int index = 0; index < (RMax-RMin+1)*(TMax-TMin+1); ++index) {
		int numberLoops = wilsonLoops[index].size();

		long_real_t result = 0.;

#pragma omp parallel for reduction(+:result)
		for (int loop = 0; loop < numberLoops; ++loop) {
			result += real(trace(wilsonLoops[index][loop].measure(environment.gaugeLinkConfiguration)));
		}

		results[index] = result/(numberLoops*numberColors);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("wilson_loops");
		for (int index = 0; index < (RMax-RMin+1)*(TMax-TMin+1); ++index) {
			std::cout << "Spatial Wilson loop (" << wilsonLoops[index][0].getR() << "," << wilsonLoops[index][0].getT() << "): " << results[index] << std::endl;
			output->push("wilson_loops");
			output->write("wilson_loops", wilsonLoops[index][0].getR());
			output->write("wilson_loops", wilsonLoops[index][0].getT());
			output->write("wilson_loops", results[index]);
			output->pop("wilson_loops");
		}
		output->pop("wilson_loops");
	}

	int numberSubSweeps = environment.configurations.get<unsigned int>("number_subsweeps_luescher");
	int sliceSize = environment.configurations.get<unsigned int>("size_slice_luescher");
	if (isOutputProcess() && (LT::glob_y % sliceSize) != 0) std::cout << "PureGaugeWilsonLoops::Warning, the Polyakov loop correlator will not work with these settings, sliceSize is not a multiple of LT::glob_y!" << std::endl;

	//Get the gauge action
	GaugeAction* action = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"), environment.configurations.get<double>("beta"));

	typedef GaugeGroup WilsonLineField[LT::glob_x][LT::glob_z][LT::glob_t][LT::glob_y/sliceSize];

	WilsonLineField *wilsonLineT = new WilsonLineField[numberSubSweeps];

	for (int numSweep = 0; numSweep < numberSubSweeps; ++numSweep) {
#pragma omp parallel for
		for (int x = 0; x < LT::glob_x; ++x) {
			for (int z = 0; z < LT::glob_z; ++z) {
				for (int t = 0; t < LT::glob_t; ++t) {
					for (int slice = 0; slice < LT::glob_y/sliceSize; ++slice) {
						int site = LT::getGlobalCoordinate(x,slice*sliceSize,z,t);
						set_to_identity(wilsonLineT[numSweep][x][z][t][slice]);
						for (int dy = 0; dy < sliceSize; ++dy) {
							wilsonLineT[numSweep][x][z][t][slice] *= environment.gaugeLinkConfiguration[site][1];
							site = extended_gauge_lattice_t::sup(site,1);
						}
					}
				}
			}
		}
		
		//Now we update the slices
		this->updateSlices(environment, action, sliceSize);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("polyakov_loop_correlator");
	}
	
	for (int R = RMin; R <= RMax; ++R) {
		long_real_t result = 0.;
		//Now we close the wilson loop
#pragma omp parallel for reduction(+:result)
		for (int x = 0; x < LT::glob_x; ++x) {
			for (int z = 0; z < LT::glob_z; ++z) {
				for (int t = 0; t < LT::glob_t; ++t) {
					TwoLinkOperator twoLinkOperator[LT::glob_y/sliceSize];
					for (int slice = 0; slice < LT::glob_y/sliceSize; ++slice) {
						twoLinkOperator[slice].setToZero();
						for (int numSweep = 0; numSweep < numberSubSweeps; ++numSweep) {
							twoLinkOperator[slice].add(wilsonLineT[numSweep][x][z][t][slice],wilsonLineT[numSweep][(x+R)%LT::glob_x][z][t][slice]);
						}
					}

					for (int slice = 0; slice < LT::glob_y/sliceSize; ++slice) {
						twoLinkOperator[slice].normalize(numberSubSweeps);
					}
					
					TwoLinkOperator twoLinkResults = twoLinkOperator[0];
					for (int slice = 1; slice < LT::glob_y/sliceSize; ++slice) {
						twoLinkResults = twoLinkResults*twoLinkOperator[slice];
					}

					for (int a = 0; a < numberColors; ++a) {
						for (int b = 0; b < numberColors; ++b) {
							result += real(twoLinkResults.at(a,a,b,b));
						}
					}
				}
			}
		}

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			output->write("polyakov_loop_correlator", result/(numberColors*numberColors*LT::glob_x*LT::glob_z*LT::glob_t));
			std::cout << "PolyakovCorrelator::Value at distance " << R << ": " << result/(numberColors*numberColors*LT::glob_x*LT::glob_z*LT::glob_t) << std::endl;
		}
		
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->pop("polyakov_loop_correlator");
	}

	delete[] wilsonLineT;
	delete action;

	environment.synchronize();
#endif
}

void PureGaugeWilsonLoops::updateSlices(environment_t& environment, GaugeAction* action, int sliceSize) {
	typedef extended_gauge_lattice_t::Layout LT;
#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif


	//Now we update, first pure gauge heatbath
#ifdef MULTITHREADING
	for (unsigned int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
#endif
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
					if ((LT::globalIndexY(site) % sliceSize) != (sliceSize - 1) && (LT::globalIndexY(site) % sliceSize) != 0) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, mu, action, action->getBeta());
#ifdef MULTITHREADING
				}
#endif
			}
#ifdef MULTITHREADING
			if (checkerboard->getColor(site,1) == color) {
#endif
				if ((LT::globalIndexY(site) % sliceSize) == 0) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, 1, action, action->getBeta());
#ifdef MULTITHREADING
			}
#endif
		}
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif

	/*//Overrelaxation
	for (int i = 0; i < 2; ++i) {
#ifdef MULTITHREADING
		for (unsigned int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
#endif
			for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
					if (checkerboard->getColor(site,mu) == color) {
#endif
						if ((LT::globalIndexY(site) % sliceSize) != (sliceSize - 1) && (LT::globalIndexY(site) % sliceSize) != 0) pureGaugeOverrelaxation->updateLink(environment.gaugeLinkConfiguration, site, mu, action);
#ifdef MULTITHREADING
					}
#endif
				}
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,1) == color) {
#endif
					if ((LT::globalIndexY(site) % sliceSize) == 0) pureGaugeOverrelaxation->updateLink(environment.gaugeLinkConfiguration, site, 1, action);
#ifdef MULTITHREADING
				}
#endif
			}
			environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
		}
#endif
	}*/
}

} /* namespace Update */
