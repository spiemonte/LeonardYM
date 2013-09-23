/*
 * PureGaugeWilsonLoops.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: spiem_01
 */

#include "PureGaugeWilsonLoops.h"
#include "GlobalOutput.h"
#include "PureGaugeUpdater.h"
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

/*//Class for computing average of links operator
class TwoLink {
	TwoLink(const TwoLink& snd) {
		for (unsigned int i = 0; i < numberColors; ++i) {
			for (unsigned int j = 0; j < numberColors; ++j) {
				for (unsigned int k = 0; k < numberColors; ++k) {
					for (unsigned int l = 0; l < numberColors; ++l) linkProduct[i][j][k][l] = snd.linkProduct[i][j][k][l];
				}
			}
		}
	}

	TwoLink& operator=(const TwoLink& snd) {
		for (unsigned int i = 0; i < numberColors; ++i) {
			for (unsigned int j = 0; j < numberColors; ++j) {
				for (unsigned int k = 0; k < numberColors; ++k) {
					for (unsigned int l = 0; l < numberColors; ++l) linkProduct[i][j][k][l] = snd.linkProduct[i][j][k][l];
				}
			}
		}
		return *this;
	}

	TwoLink operator*(const TwoLink& snd) {
		TwoLink ris;
		for (unsigned int alpha = 0; alpha < numberColors; ++alpha) {
			for (unsigned int beta = 0; beta < numberColors; ++beta) {
				for (unsigned int gamma = 0; gamma < numberColors; ++gamma) {
					for (unsigned int delta = 0; delta < numberColors; ++delta) {
						ris.linkProduct[alpha][beta][gamma][delta] = 0.;
						for (unsigned int lambda = 0; lambda < numberColors; ++lambda) {
							for (unsigned int epsilon = 0; epsilon < numberColors; ++epsilon) {
								ris.linkProduct[alpha][beta][gamma][delta] += linkProduct[alpha][lambda][gamma][epsilon]*snd.linkProduct[lambda][beta][epsilon][delta];
							}
						}
					}
				}
			}
		}
		return ris;
	}

	std::complex<real_t>& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
		return linkProduct[i][j][k][l];
	}

private:
	std::complex<real_t> linkProduct[numberColors][numberColors][numberColors][numberColors];
};

class SpatialTwoLinkStorage {
public:
	SpatialTwoLinkStorage(int _z, int _t, int _x1, int _x2, int _y1, int _y2) :
	 z(_z), t(_t), x1(_x1), x2(_x2), y1(_y1), y2(_y2) { }

	TwoLink getTwoLink() {
		TwoLink result;
		for (unsigned int alpha = 0; alpha < numberColors; ++alpha) {
			for (unsigned int beta = 0; beta < numberColors; ++beta) {
				for (unsigned int gamma = 0; gamma < numberColors; ++gamma) {
					for (unsigned int delta = 0; delta < numberColors; ++delta) {
						std::complex<long_real_t> mean(0.,0.);
						for (int i = 0; i < links1.size(); ++i) {
							mean += conj(links1[i].at(alpha,beta))*(links2[i].at(gamma,delta));
						}
						result(alpha,beta,gamma,delta) = std::complex<real_t>(real(mean)/links1.size(),imag(mean)/links1.size());
					}
				}
			}
		}
		return result;
	}

	void reset() {
		links1.clear();
		links2.clear();
	}

	void addMeasure(const fundamental_lattice_t& lattice) {
		GaugeGroup link;
		set_to_identity(link);
		typedef fundamental_lattice_t::Layout Layout;
		for (unsigned int y = y1; y < y2; ++y) {
			int site = Layout::index(x1,y,z,t);
			link = link*lattice[site][1];
		}
		links1.push_back(link);
		set_to_identity(link);
		for (unsigned int y = y1; y < y2; ++y) {
			int site = Layout::index(x2,y,z,t);
			link = link*lattice[site][1];
		}
		links2.push_back(link);
	}

private:
	int z;
	int t;
	int x1;
	int x2;
	int y1;
	int y2;

	std::vector<GaugeGroup> links1;
	std::vector<GaugeGroup> links2;
};

class ActiveLinks {
public:
	ActiveLinks() { }

	void addBand(int y1, int y2) {

	}
private:
	bool (*activeLinks)[4];
};*/

PureGaugeWilsonLoops::PureGaugeWilsonLoops() : wilsonLoops(0), results(0) { }

PureGaugeWilsonLoops::~PureGaugeWilsonLoops() {
	if (results) delete[] results;
	if (wilsonLoops) delete[] wilsonLoops;
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


#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

	real_t beta = environment.configurations.get<real_t>("beta");

	int numberSubSweeps = environment.configurations.get<unsigned int>("number_subsweeps_luescher");
	int sliceSize = environment.configurations.get<unsigned int>("size_slice_luescher");
	if (isOutputProcess() && (LT::glob_y % sliceSize) != 0) std::cout << "PureGaugeWilsonLoops::Warning, the Polyakov loop correlator will not work with these settings, sliceSize is not a multiple of LT::glob_y!" << std::endl;

	//Get the gauge action
	GaugeAction* action = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

	PureGaugeUpdater* pureGaugeUpdater = new PureGaugeUpdater();

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
		
		//Now we update
#ifdef MULTITHREADING
		for (unsigned int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
#endif
			for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
					if (checkerboard->getColor(site,mu) == color) {
#endif
						if ((LT::globalIndexY(site) % sliceSize) != (sliceSize - 1) && (LT::globalIndexY(site) % sliceSize) != 0) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta);
#ifdef MULTITHREADING
					}
#endif
				}
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,1) == color) {
#endif
					if ((LT::globalIndexY(site) % sliceSize) == 0) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, 1, action, beta);
#ifdef MULTITHREADING
				}
#endif
			}
			environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
		}
#endif
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
	delete pureGaugeUpdater;
	delete action;















	/*GaugeGroup ****wilsonLineT = new GaugeGroup***[numberSubSweeps];
	for (int i = 0; i < numberSubSweeps; ++i) {
		wilsonLineT[i] = new GaugeGroup**[LT::glob_x];
		for (int x = 0; x < LT::glob_x; ++x) {
			wilsonLineT[i][x] = new GaugeGroup*[LT::glob_z];
			for (int z = 0; z < LT::glob_z; ++z) {
				wilsonLineT[i][x][z] = new GaugeGroup[LT::glob_t];
			}
		}
	}

	for (unsigned int T = TMin; T <= TMax; ++T) {
		
		for (int numSweep = 0; numSweep < numberSubSweeps; ++numSweep) {
#pragma omp parallel for
			for (int x = 0; x < LT::glob_x; ++x) {
				for (int z = 0; z < LT::glob_z; ++z) {
					for (int t = 0; t < LT::glob_t; ++t) {
						int site = LT::getGlobalCoordinate(x,0,z,t);
						set_to_identity(wilsonLineT[numSweep][x][z][t]);
						for (int dy = 0; dy < T; ++dy) {
							wilsonLineT[numSweep][x][z][t] *= environment.gaugeLinkConfiguration[site][1];
							site = extended_gauge_lattice_t::sup(site,1);
						}
					}
				}
			}
			//Now we update
#ifdef MULTITHREADING
			for (unsigned int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
#endif
				for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
						if (checkerboard->getColor(site,mu) == color) {
#endif
							if (LT::globalIndexY(site) > 0 && LT::globalIndexY(site) < T) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta);
#ifdef MULTITHREADING
						}
#endif
					}
#ifdef MULTITHREADING
					if (checkerboard->getColor(site,1) == color) {
#endif
						if (LT::globalIndexY(site) == 0) pureGaugeUpdater->updateLink(environment.gaugeLinkConfiguration, site, 1, action, beta);
#ifdef MULTITHREADING
					}
#endif
				}
				environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
			}
#endif
		}

		for (unsigned int R = RMin; R <= RMax; ++R) {
			long_real_t result = 0.;
			//Now we close the wilson loop
#pragma omp parallel for reduction(+:result)
			for (int x = 0; x < LT::glob_x; ++x) {
				for (int z = 0; z < LT::glob_z; ++z) {
					for (int t = 0; t < LT::glob_t; ++t) {
						TwoLinkOperator twoLinkOperator;
						twoLinkOperator.setToZero();
						for (int numSweep = 0; numSweep < numberSubSweeps; ++numSweep) {
							twoLinkOperator.add(wilsonLineT[numSweep][x][z][t],wilsonLineT[numSweep][(x+R)%LT::glob_x][z][t]);
						}
						
						GaugeGroup U1;
						set_to_identity(U1);
						int site = LT::getGlobalCoordinate(x,0,z,t);
						for (int dx = 0; dx < R; ++dx) {
							U1 *= environment.gaugeLinkConfiguration[site][0];
							site = extended_gauge_lattice_t::sup(site,0);
						}
						GaugeGroup U2;
						set_to_identity(U2);
						site = LT::getGlobalCoordinate(x,T,z,t);
						for (int dx = 0; dx < R; ++dx) {
							U2 *= environment.gaugeLinkConfiguration[site][0];
							site = extended_gauge_lattice_t::sup(site,0);
						}
						GaugeGroup hU1 = htrans(U1);
						for (int a = 0; a < numberColors; ++a) {
							for (int b = 0; b < numberColors; ++b) {
								for (int c = 0; c < numberColors; ++c) {
									for (int d = 0; d < numberColors; ++d) {
										result += real(twoLinkOperator.at(a,b,c,d)*U2.at(b,c)*hU1.at(d,a))/numberSubSweeps;
									}
								}
							}
						}
					}
				}
			}
			wilsonLoopResult.push_back(result/(LT::glob_x*LT::glob_z*LT::glob_t*numberColors));
			wilsonLoopR.push_back(R);
			wilsonLoopT.push_back(T);
		}

		//We wrap the lattice to reduce the cross correlations
		//since the periodic boundary conditions, nothing should change

		//We do a swap
		for(int site = 0; site < (LT::localsize); ++site){
			for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = environment.gaugeLinkConfiguration[site][mu];
		}
		swaplinkconfig.updateHalo();
		//We wrap
		for(int site = 0; site < (LT::localsize); ++site){
			for (unsigned int mu = 0; mu < 4; ++mu) environment.gaugeLinkConfiguration[site][mu] = swaplinkconfig[extended_gauge_lattice_t::sup(site,1)][mu];
		}
		environment.gaugeLinkConfiguration.updateHalo();
	}

	for (int i = 0; i < numberSubSweeps; ++i) {
		for (int x = 0; x < LT::glob_x; ++x) {
			for (int z = 0; z < LT::glob_z; ++z) {
				delete[] wilsonLineT[i][x][z];
			}
			delete[] wilsonLineT[i][x];
		}
		delete[] wilsonLineT[i];
	}
	delete[] wilsonLineT;*/

	environment.synchronize();
#endif
}

} /* namespace Update */
