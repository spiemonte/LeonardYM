/*
 * WilsonFlow.cpp
 *
 *  Created on: Oct 23, 2013
 *      Author: spiem_01
 */

#include "WilsonFlow.h"
#include "ImprovedGaugeAction.h"
#include "WilsonGaugeAction.h"
#include "GlobalOutput.h"
#include "Plaquette.h"

namespace Update {

WilsonFlow::WilsonFlow() { }

WilsonFlow::~WilsonFlow() { }

void WilsonFlow::execute(environment_t& environment) {
	extended_gauge_lattice_t initialLattice = environment.gaugeLinkConfiguration, finalLattice;
	std::string flow_type = environment.configurations.get<std::string>("flow_type");
	GaugeAction* action;
	if (flow_type == "Wilson") {
		action = new WilsonGaugeAction(2*numberColors);
	}
	else if (flow_type == "Symanzik") {
		action = new ImprovedGaugeAction(2*numberColors, 1);
	}
	else {
		std::cout << "WilsonFlow::Flow type " << flow_type << " unknown! (allowed types: Wilson/Symanzik)" << std::endl;
		exit(1);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("wilson_flow");
	}

	real_t step = environment.configurations.get<real_t>("flow_step");
	real_t flow_time = environment.configurations.get<real_t>("flow_time");
	unsigned int integration_intervals = environment.configurations.get<unsigned int>("flow_integration_intervals");

	for (real_t t = 0; t < flow_time; t += step) {
		real_t energy = this->measureEnergy(initialLattice);
		std::cout << "WilsonFlow::t*t*Energy at t " << t << ": " << t*t*energy << std::endl;
		this->integrate(initialLattice, finalLattice, action, step, integration_intervals);
		initialLattice = finalLattice;

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			output->push("wilson_flow");
			output->write("wilson_flow", t);
			output->write("wilson_flow", t*t*energy);
			output->pop("wilson_flow");
		}
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->pop("wilson_flow");
	}

	delete action;
}

void WilsonFlow::integrate(const extended_gauge_lattice_t& initialLattice, extended_gauge_lattice_t& finalLattice, GaugeAction* action, real_t time, int nSteps) {
	/* notations from hep-lat/1203.4469 */
	real_t dt = time/nSteps;
	extended_gauge_lattice_t Z0, Z1, Z2;
	extended_gauge_lattice_t X0 = initialLattice, X1, X2;
	for (int step = 0; step < nSteps; ++step) {
		this->getForce(X0, Z0, action);

#pragma omp parallel for
		for (int site = 0; site < Z0.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				X1[site][mu] = this->exponential(X0[site][mu], (1/4.)*Z0[site][mu], dt);
			}
		}
		X1.updateHalo();

		this->getForce(X1, Z1, action);

#pragma omp parallel for
		for (int site = 0; site < Z1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				X2[site][mu] = this->exponential(X1[site][mu], (8./9.)*Z1[site][mu] - (17./36.)*Z0[site][mu], dt);
			}
		}
		X2.updateHalo();

		this->getForce(X2, Z2, action);

#pragma omp parallel for
		for (int site = 0; site < Z1.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				X0[site][mu] = this->exponential(X2[site][mu], (3./4.)*Z2[site][mu] - (8./9.)*Z1[site][mu] + (17./36.)*Z0[site][mu], dt);
			}
		}
		X0.updateHalo();
	}

	finalLattice = X0;
}

void WilsonFlow::getForce(const extended_gauge_lattice_t& lattice, extended_gauge_lattice_t& force, GaugeAction* action) {
#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			force[site][mu] = -action->force(lattice, site, mu);
		}
	}
}

GaugeGroup WilsonFlow::exponential(const GaugeGroup& link, const GaugeGroup& force, real_t epsilon) {
#ifdef EIGEN
	Eigen::ComplexEigenSolver<GaugeGroup> es(force);
	GaugeGroup update;
	update.zeros();
	for (int i = 0; i < numberColors; ++i) {
		update.at(i,i) = exp(epsilon*es.eigenvalues()[i]);
	}
	GaugeGroup m = es.eigenvectors();
	GaugeGroup updatenew = m * update * htrans(m);
#endif
#ifdef ARMADILLO
	FundamentalVector eigval;
	GaugeGroup eigvec;
	arma::eig_gen(eigval, eigvec, force);
	GaugeGroup update;
	set_to_zero(update);
	for (unsigned int i = 0; i < numberColors; ++i) {
		update.at(i,i) = exp(epsilon*eigval[i]);
	}
	GaugeGroup updatenew = eigvec * update * htrans(eigvec);
#endif
#ifdef MATRIXTOOLKIT
	FundamentalVector eigval;
	GaugeGroup eigvec;
	matrix_toolkit::eigensystem(eigval, eigvec, force);
	GaugeGroup update;
	set_to_zero(update);
	for (unsigned int i = 0; i < numberColors; ++i) {
		update.at(i,i) = exp(epsilon*eigval[i]);
	}
	GaugeGroup updatenew = eigvec * update * htrans(eigvec);
#endif
	return updatenew*link;
}

real_t WilsonFlow::measureEnergy(const extended_gauge_lattice_t& _lattice) {
	typedef extended_fermion_lattice_t LT;
	typedef extended_fermion_lattice_t::Layout Layout;
	long_real_t energy = 0.;
#pragma omp parallel for reduction(+:energy)
	for (int site = 0; site < _lattice.localsize; ++site) {
		GaugeGroup tmpF[6];
		tmpF[0] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 1)][1])*(_lattice[LT::sdn(LT::sdn(site, 0), 1)][0])*(_lattice[LT::sdn(site, 1)][1]) + htrans(_lattice[LT::sdn(site, 1)][1])*(_lattice[LT::sdn(site, 1)][0])*(_lattice[LT::sup(LT::sdn(site, 1), 0)][1])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][1])*htrans(_lattice[LT::sup(site, 1)][0])*htrans(_lattice[site][1]) + (_lattice[site][1])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 1)][0])*htrans(_lattice[LT::sdn(site, 0)][1])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[1] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 2)][2])*(_lattice[LT::sdn(LT::sdn(site, 0), 2)][0])*(_lattice[LT::sdn(site, 2)][2]) + htrans(_lattice[LT::sdn(site, 2)][2])*(_lattice[LT::sdn(site, 2)][0])*(_lattice[LT::sup(LT::sdn(site, 2), 0)][2])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][2])*htrans(_lattice[LT::sup(site, 2)][0])*htrans(_lattice[site][2]) + (_lattice[site][2])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 2)][0])*htrans(_lattice[LT::sdn(site, 0)][2])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[2] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 0), 3)][0])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][0])*(_lattice[LT::sup(LT::sdn(site, 3), 0)][3])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][3])*htrans(_lattice[LT::sup(site, 3)][0])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 3)][0])*htrans(_lattice[LT::sdn(site, 0)][3])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[3] = htrans(_lattice[LT::sdn(site, 1)][1])*htrans(_lattice[LT::sdn(LT::sdn(site, 1), 2)][2])*(_lattice[LT::sdn(LT::sdn(site, 1), 2)][1])*(_lattice[LT::sdn(site, 2)][2]) + htrans(_lattice[LT::sdn(site, 2)][2])*(_lattice[LT::sdn(site, 2)][1])*(_lattice[LT::sup(LT::sdn(site, 2), 1)][2])*htrans(_lattice[site][1]) + (_lattice[site][1])*(_lattice[LT::sup(site, 1)][2])*htrans(_lattice[LT::sup(site, 2)][1])*htrans(_lattice[site][2]) + (_lattice[site][2])*htrans(_lattice[LT::sup(LT::sdn(site, 1), 2)][1])*htrans(_lattice[LT::sdn(site, 1)][2])*(_lattice[LT::sdn(site, 1)][1]);
		tmpF[4] = htrans(_lattice[LT::sdn(site, 1)][1])*htrans(_lattice[LT::sdn(LT::sdn(site, 1), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 1), 3)][1])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][1])*(_lattice[LT::sup(LT::sdn(site, 3), 1)][3])*htrans(_lattice[site][1]) + (_lattice[site][1])*(_lattice[LT::sup(site, 1)][3])*htrans(_lattice[LT::sup(site, 3)][1])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 1), 3)][1])*htrans(_lattice[LT::sdn(site, 1)][3])*(_lattice[LT::sdn(site, 1)][1]);
		tmpF[5] = htrans(_lattice[LT::sdn(site, 2)][2])*htrans(_lattice[LT::sdn(LT::sdn(site, 2), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 2), 3)][2])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][2])*(_lattice[LT::sup(LT::sdn(site, 3), 2)][3])*htrans(_lattice[site][2]) + (_lattice[site][2])*(_lattice[LT::sup(site, 2)][3])*htrans(_lattice[LT::sup(site, 3)][2])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 2), 3)][2])*htrans(_lattice[LT::sdn(site, 2)][3])*(_lattice[LT::sdn(site, 2)][2]);
		for (unsigned int i = 0; i < 6; ++i) {
			//Manual antialiasing, error of eigen!
			GaugeGroup antialias = tmpF[i];
			tmpF[i] = (1./8.)*(antialias - htrans(antialias));
			energy += 4.*(1/4.)*real(trace(tmpF[i]*tmpF[i]));
		}
	}
	reduceAllSum(energy);
	return -energy/Layout::globalVolume;
}

} /* namespace Update */
