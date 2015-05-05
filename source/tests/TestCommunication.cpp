/*
 * TestCommunication.cpp
 *
 *  Created on: Feb 25, 2013
 *      Author: spiem_01
 */

#include "TestCommunication.h"
#include "AlgebraUtils.h"

namespace Update {

TestCommunication::TestCommunication() { }

TestCommunication::~TestCommunication() { }

void TestCommunication::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t LT;
	extended_gauge_lattice_t lattice = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t swaplinkconfig;
	double energyWilson1 = 0.;

	//First we measure the standard wilson energy
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += lattice[LT::sup(site,mu)][nu]*adj(lattice[LT::sup(site,nu)][mu])*adj(lattice[site][nu]);
			}
			energyWilson1 += -real(trace(lattice[site][mu]*plaqs));
		}
	}
	reduceAllSum(energyWilson1);

	//Then we measure the improved gauge action before
	double energyImproved1 = 0.;

	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += (5./3.)*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
				plaqs -= (1./(12.))*(lattice[LT::sup(site, mu)][mu])*(lattice[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(lattice[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[site][nu]);
				plaqs -= (1./(12.))*(lattice[LT::sup(site, mu)][nu])*(lattice[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(lattice[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][nu])*htrans(lattice[site][nu]);
			}
			energyImproved1 += -real(trace(lattice[site][mu]*plaqs));
		}
	}
	reduceAllSum(energyImproved1);

	//Then we measure the "staple improved energy", the energy given by the staple multiplied with a link.
	//It contains more only the "negative" rectangles and squares
	double stapleEnergy1 = 0.;
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup ris;
			set_to_zero(ris);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					ris += (5./3.)*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
					ris += (5./3.)*htrans(lattice[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(lattice[LT::sdn(site,nu)][mu])*lattice[LT::sdn(site,nu)][nu];
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][mu])*(lattice[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(lattice[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[site][nu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][nu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[LT::sup(LT::sdn(site, mu), nu)][mu])*htrans(lattice[LT::sdn(site, mu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][mu])*htrans(lattice[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu])*htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][mu])*htrans(lattice[LT::sdn(site, nu)][mu])*(lattice[LT::sdn(site, nu)][nu]);
					ris -= (1./(12.))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sdn(site, nu)][mu])*htrans(lattice[LT::sdn(LT::sdn(site, mu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, mu), nu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][nu])*(lattice[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(lattice[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][nu])*htrans(lattice[site][nu]);
					ris -= (1./(12.))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu])*htrans(lattice[LT::sdn(LT::sdn(site, nu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, nu), nu)][nu])*(lattice[LT::sdn(site, nu)][nu]);
				}
			}
			stapleEnergy1 += -real(trace(lattice[site][mu]*ris));
		}
	}
	reduceAllSum(stapleEnergy1);

	//Now we simply wrap around of one site in the temporal direction the linkconfig,
	//since the periodic boundary conditions, nothing should change

	//We do a swap
	for(int site = 0; site < (lattice.localsize); ++site){
		for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = lattice[site][mu];
	}
	swaplinkconfig.updateHalo();
	//We wrap
	for(int site = 0; site < (lattice.localsize); ++site){
		for (unsigned int mu = 0; mu < 4; ++mu) lattice[site][mu] = swaplinkconfig[LT::sup(site,0)][mu];
	}
	lattice.updateHalo();

	//Then we measure the standard gauge action after
	double energyWilson2 = 0.;
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += lattice[LT::sup(site,mu)][nu]*adj(lattice[LT::sup(site,nu)][mu])*adj(lattice[site][nu]);
			}
			energyWilson2 += -real(trace(lattice[site][mu]*plaqs));
		}
	}
	reduceAllSum(energyWilson2);

	if (isOutputProcess()) std::cout << "TestCommunication::Difference in Wilson gauge energy: " << energyWilson2 - energyWilson1 << std::endl;

	//Then we measure the improved gauge action after
	double energyImproved2 = 0.;

	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup plaqs;
			set_to_zero(plaqs);
			for (unsigned int nu = mu + 1; nu < 4; ++nu) {
				plaqs += (5./3.)*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
				plaqs -= (1./(12.))*(lattice[LT::sup(site, mu)][mu])*(lattice[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(lattice[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[site][nu]);
				plaqs -= (1./(12.))*(lattice[LT::sup(site, mu)][nu])*(lattice[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(lattice[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][nu])*htrans(lattice[site][nu]);
			}
			energyImproved2 += -real(trace(lattice[site][mu]*plaqs));
		}
	}
	reduceAllSum(energyImproved2);

	if (isOutputProcess()) std::cout << "TestCommunication::Difference in Improved gauge energy: " << energyImproved2 - energyImproved1 << std::endl;

	double stapleEnergy2 = 0.;
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup ris;
			set_to_zero(ris);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					ris += (5./3.)*lattice[LT::sup(site,mu)][nu]*htrans(lattice[LT::sup(site,nu)][mu])*htrans(lattice[site][nu]);
					ris += (5./3.)*htrans(lattice[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(lattice[LT::sdn(site,nu)][mu])*lattice[LT::sdn(site,nu)][nu];

					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][mu])*(lattice[LT::sup(LT::sup(site, mu), mu)][nu])*htrans(lattice[LT::sup(LT::sup(site, mu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[site][nu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][nu])*htrans(lattice[LT::sup(site, nu)][mu])*htrans(lattice[LT::sup(LT::sdn(site, mu), nu)][mu])*htrans(lattice[LT::sdn(site, mu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][mu])*htrans(lattice[LT::sdn(LT::sup(LT::sup(site, mu), mu), nu)][nu])*htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][mu])*htrans(lattice[LT::sdn(site, nu)][mu])*(lattice[LT::sdn(site, nu)][nu]);
					ris -= (1./(12.))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sdn(site, nu)][mu])*htrans(lattice[LT::sdn(LT::sdn(site, mu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, mu), nu)][nu])*(lattice[LT::sdn(site, mu)][mu]);
					ris -= (1./(12.))* (lattice[LT::sup(site, mu)][nu])*(lattice[LT::sup(LT::sup(site, mu), nu)][nu])*htrans(lattice[LT::sup(LT::sup(site, nu), nu)][mu])*htrans(lattice[LT::sup(site, nu)][nu])*htrans(lattice[site][nu]);
					ris -= (1./(12.))* htrans(lattice[LT::sup(LT::sdn(site, nu), mu)][nu])*htrans(lattice[LT::sup(LT::sdn(LT::sdn(site, nu), nu), mu)][nu])*htrans(lattice[LT::sdn(LT::sdn(site, nu), nu)][mu])*(lattice[LT::sdn(LT::sdn(site, nu), nu)][nu])*(lattice[LT::sdn(site, nu)][nu]);
				}
			}
			stapleEnergy2 += -real(trace(lattice[site][mu]*ris));
		}
	}
	reduceAllSum(stapleEnergy2);

	if (isOutputProcess()) std::cout <<"TestCommunication::Difference \"improved staple energy\": "<< stapleEnergy2 - stapleEnergy1 << std::endl;

	typedef reduced_dirac_vector_t DV;
	reduced_dirac_vector_t test1, test2, test3;
	AlgebraUtils::generateRandomVector(test1);
	long_real_t before = AlgebraUtils::squaredNorm(test1);

	for (int site = 0; site < test1.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			test2[site][mu] = test1[DV::sup(site,2)][mu];
		}
	}

	long_real_t after1 = AlgebraUtils::squaredNorm(test2);

	for (int site = 0; site < test1.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			test3[site][mu] = test1[DV::sdn(site,2)][mu];
		}
	}

	long_real_t after2 = AlgebraUtils::squaredNorm(test3);

	if (isOutputProcess()) std::cout << "TestCommunication::Test communication dirac vector: " << before - after1 << " " << before - after2 << std::endl;
}

} /* namespace Update */
