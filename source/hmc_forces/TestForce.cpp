/*
 * TestForce.cpp
 *
 *  Created on: Jun 4, 2012
 *      Author: spiem_01
 */

#include "TestForce.h"
#include "ToString.h"
#include "Plaquette.h"

namespace Update {

TestForce::TestForce() { }

TestForce::~TestForce() { }

void TestForce::testForce(const environment_t& env, Energy* energy, const GaugeGroup& force, int site, int mu) {
#if NUMCOLORS == 2
	//First create a "copy" environment
	tmp_env = env;

	//then evaluate the force
	real_t epsilon = 0.01;
	while (epsilon > 0.00001) {
		//The raw force
		matrix2x2_t raw_force;
		raw_force(0,0) = 0.;
		raw_force(0,1) = 0.;
		raw_force(1,0) = 0.;
		raw_force(1,1) = 0.;

		//First component
		matrix2x2_t expm;
		expm(0,0) = complex(cos(epsilon), 0.);
		expm(0,1) = complex(0., sin(epsilon));
		expm(1,0) = complex(0., sin(epsilon));
		expm(1,1) = complex(cos(epsilon), 0.);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*(tmp_env.gaugeLinkConfiguration[site][mu]);
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		long_real_t component0plus = energy->energy(tmp_env);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		expm(0,0) = complex(cos(epsilon), 0.);
		expm(0,1) = complex(0., -sin(epsilon));
		expm(1,0) = complex(0., -sin(epsilon));
		expm(1,1) = complex(cos(epsilon), 0.);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		long_real_t component0minus = energy->energy(tmp_env);
		raw_force(0,0) += 0.;
		raw_force(0,1) += complex(0., (component0plus-component0minus)/(8*epsilon));
		raw_force(1,0) += complex(0., (component0plus-component0minus)/(8*epsilon));
		raw_force(1,1) += 0.;
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		if (isOutputProcess()) std::cout << "Force:: This must be the sensibility: " << (component0plus-component0minus) << std::endl;

		//Second component
		expm(0,0) = complex(cos(epsilon), 0.);
		expm(0,1) = complex(sin(epsilon), 0.);
		expm(1,0) = complex(-sin(epsilon), 0.);
		expm(1,1) = complex(cos(epsilon), 0.);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		component0plus = energy->energy(tmp_env);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		expm(0,0) = complex(cos(epsilon), 0.);
		expm(0,1) = complex(-sin(epsilon), 0.);
		expm(1,0) = complex(sin(epsilon), 0.);
		expm(1,1) = complex(cos(epsilon), 0.);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		component0minus = energy->energy(tmp_env);
		raw_force(0,0) += 0.;
		raw_force(0,1) += complex(+(component0plus-component0minus)/(8*epsilon), 0.);
		raw_force(1,0) += complex(-(component0plus-component0minus)/(8*epsilon), 0.);
		raw_force(1,1) += 0.;
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();

		if (isOutputProcess()) std::cout << "Force:: This must be the sensibility: " << (component0plus-component0minus) << std::endl;

		//Third component
		expm(0,0) = complex(cos(epsilon), sin(epsilon));
		expm(0,1) = complex(0., 0.);
		expm(1,0) = complex(0., 0.);
		expm(1,1) = complex(cos(epsilon), -sin(epsilon));
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		component0plus = energy->energy(tmp_env);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		expm(0,0) = complex(cos(epsilon), -sin(epsilon));
		expm(0,1) = complex(0., 0.);
		expm(1,0) = complex(0., 0.);
		expm(1,1) = complex(cos(epsilon), sin(epsilon));
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		component0minus = energy->energy(tmp_env);
		raw_force(0,0) += complex(0., +(component0plus-component0minus)/(8*epsilon));
		raw_force(0,1) += 0.;
		raw_force(1,0) += 0.;
		raw_force(1,1) += complex(0., -(component0plus-component0minus)/(8*epsilon));
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();

		if (isOutputProcess()) std::cout << "Force:: This must be the sensibility: " << (component0plus-component0minus) << std::endl;

		//Now we can test
		if (isOutputProcess()) {
			std::cout << "Test of the force (numerical - provided), epsilon: " << epsilon << std::endl;
			std::cout << "Test of the force, component (0,0): " << raw_force(0,0) << " " << force(0,0) << std::endl;
			std::cout << "Test of the force, component (0,1): " << raw_force(0,1) << " " << force(0,1) << std::endl;
			std::cout << "Test of the force, component (1,0): " << raw_force(1,0) << " " << force(1,0) << std::endl;
			std::cout << "Test of the force, component (1,1): " << raw_force(1,1) << " " << force(1,1) << std::endl;
		}
		epsilon = epsilon/10.;
	}
#endif
}

void TestForce::genericTestForce(const environment_t& env, Energy* energy, const GaugeGroup& force, int site, int mu) {
	typedef extended_gauge_lattice_t LT;
	//First create a "copy" environment
	tmp_env = env;

	//then evaluate the force
	real_t epsilon = 0.01;
	while (epsilon > 0.00001) {
		//First component
		GaugeGroup expm;
		set_to_zero(expm);

		//Third component
		for (unsigned int i = 0; i < numberColors - 1; ++i) {
			expm(i,i) = complex(cos(sqrt(2./static_cast<double>(numberColors*(numberColors-1)))*epsilon), sin(sqrt(2./static_cast<double>(numberColors*(numberColors-1)))*epsilon));
		}
		expm(numberColors-1,numberColors-1) = complex(cos(sqrt((2.*(numberColors - 1))/static_cast<double>(numberColors))*epsilon), -sin(sqrt((2.*(numberColors - 1))/static_cast<double>(numberColors))*epsilon));

		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();
		long_real_t component0plus = energy->energy(tmp_env);
		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();

		for (unsigned int i = 0; i < numberColors - 1; ++i) {
			expm(i,i) = complex(cos(sqrt(2./static_cast<double>(numberColors*(numberColors-1)))*epsilon), -sin(sqrt(2./static_cast<double>(numberColors*(numberColors-1)))*epsilon));
		}
		expm(numberColors-1,numberColors-1) = complex(cos(sqrt((2.*(numberColors - 1))/static_cast<double>(numberColors))*epsilon), sin(sqrt((2.*(numberColors - 1))/static_cast<double>(numberColors))*epsilon));

		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = expm*tmp_env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();

		long_real_t component0minus = energy->energy(tmp_env);

		if (isOutputProcess()) tmp_env.gaugeLinkConfiguration[site][mu] = env.gaugeLinkConfiguration[site][mu];
		tmp_env.gaugeLinkConfiguration.updateHalo();
		tmp_env.synchronize();

		if (isOutputProcess()) std::cout << "Force:: This must be the sensibility: " << (component0plus-component0minus) << std::endl;

		//Now we can test
		if (isOutputProcess()) {
			std::cout << "Test of the force (numerical - provided), epsilon: " << epsilon << std::endl;
			std::cout << "Test of the force, component (nc-1,nc-1): " << -sqrt((2.*(numberColors - 1))/static_cast<double>(numberColors))*(component0plus-component0minus)/(8*epsilon) << " " << imag(force(numberColors-1,numberColors-1)) << std::endl;
		}
		epsilon = epsilon/10.;
	}
}

} /* namespace Update */
