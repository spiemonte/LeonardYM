/*
 * PureGaugeUpdater.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: spiem_01
 */

#include "PureGaugeUpdater.h"
#include "MatrixTypedef.h"
#include "Checkerboard.h"
#include "RandomSeed.h"
#ifndef PI
#define PI 3.14159265358979323846264338327950288419
#endif
#include <omp.h>

namespace Update {

#ifdef MULTITHREADING
PureGaugeUpdater::PureGaugeUpdater() {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdater::PureGaugeUpdater() : randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif

#ifdef MULTITHREADING
PureGaugeUpdater::PureGaugeUpdater(const PureGaugeUpdater& copy) : LatticeSweep(copy) {
	randomGenerator = new random_generator_t*[omp_get_max_threads()];
	randomUniform = new random_uniform_generator_t*[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		randomGenerator[i] = new random_generator_t(RandomSeed::randomSeed());
		randomUniform[i] = new random_uniform_generator_t(RandomSeed::getRandomNumberGenerator(*randomGenerator[i]));
	}
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdater::PureGaugeUpdater(const PureGaugeUpdater& copy) : randomGenerator(RandomSeed::randomSeed()), randomUniform(RandomSeed::getRandomNumberGenerator(randomGenerator)) { }
#endif

#ifdef MULTITHREADING
PureGaugeUpdater::~PureGaugeUpdater() {
	for (int i = 0; i < omp_get_max_threads(); ++i) {
		delete randomGenerator[i];
		delete randomUniform[i];
	}
	delete[] randomGenerator;
	delete[] randomUniform;
}
#endif
#ifndef MULTITHREADING
PureGaugeUpdater::~PureGaugeUpdater() { }
#endif

void PureGaugeUpdater::execute(environment_t & environment) {
	typedef extended_gauge_lattice_t::Layout Layout;
	real_t beta = environment.configurations.get<double>("beta");

	//Get the gauge action
	GaugeAction* action = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

	
#ifdef MULTITHREADING
	Checkerboard* checkerboard = Checkerboard::getInstance();
#endif

#ifndef ENABLE_MPI
#ifdef MULTITHREADING
	for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
#endif
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
		for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
#ifdef MULTITHREADING
				if (checkerboard->getColor(site,mu) == color) {
#endif
					this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta);
#ifdef MULTITHREADING
				}
#endif
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
#ifdef MULTITHREADING
	}
#endif
#endif

	//We suppose that always MPI+MTH
#ifdef ENABLE_MPI
	for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
		for (int color = 0; color < checkerboard->getNumberLoops(); ++color) {
			if (processor == Layout::this_processor) {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
				for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta);
						}
					}
				}
			}
			else {
#pragma omp parallel for //shared(beta, color, environment) firstprivate(action, checkerboard) default(none) schedule(dynamic)
				for (int site = environment.gaugeLinkConfiguration.sharedsize; site < environment.gaugeLinkConfiguration.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						if (checkerboard->getColor(site,mu) == color) {
							this->updateLink(environment.gaugeLinkConfiguration, site, mu, action, beta);
						}
					}
				}
			}
		}
		environment.gaugeLinkConfiguration.updateHalo();
	}
#endif

	environment.synchronize();
	delete action;
}

#ifdef MULTITHREADING
real_t PureGaugeUpdater::generate_radius(real_t b) {
	real_t x1 = log((*randomUniform[omp_get_thread_num()])()), x2 = log((*randomUniform[omp_get_thread_num()])()), x3 = pow(cos(2.*PI*(*randomUniform[omp_get_thread_num()])()), 2.);
	real_t s = 1. + b*(x1+x2*x3);
	real_t r = (*randomUniform[omp_get_thread_num()])();
	while ((1.+s-2.*r*r) < 0) {
		x1 = log((*randomUniform[omp_get_thread_num()])());
		x2 = log((*randomUniform[omp_get_thread_num()])());
		x3 = pow(cos(2.*PI*(*randomUniform[omp_get_thread_num()])()), 2.);
		s = 1. + b*(x1+x2*x3);
		r = (*randomUniform[omp_get_thread_num()])();
	}
	return s;
}
#endif
#ifndef MULTITHREADING
real_t PureGaugeUpdater::generate_radius(real_t b) {
	real_t x1 = log(randomUniform()), x2 = log(randomUniform()), x3 = pow(cos(2.*PI*randomUniform()), 2.);
	real_t s = 1. + b*(x1+x2*x3);
	real_t r = randomUniform();
	while ((1.+s-2.*r*r) < 0) {
		x1 = log(randomUniform());
		x2 = log(randomUniform());
		x3 = pow(cos(2.*PI*randomUniform()), 2.);
		s = 1. + b*(x1+x2*x3);
		r = randomUniform();
	}
	return s;
}
#endif

#ifdef MULTITHREADING
void PureGaugeUpdater::generate_vector(real_t radius, real_t& u1, real_t& u2, real_t& u3) {
	real_t phi = 2.*PI*(*randomUniform[omp_get_thread_num()])();
	real_t theta = acos(2.*(*randomUniform[omp_get_thread_num()])()-1.);
	u1 = radius*sin(theta)*cos(phi);
	u2 = radius*sin(theta)*sin(phi);
	u3 = radius*cos(theta);
}
#endif
#ifndef MULTITHREADING
void PureGaugeUpdater::generate_vector(real_t radius, real_t& u1, real_t& u2, real_t& u3) {
	real_t phi = 2.*PI*randomUniform();
	real_t theta = acos(2.*randomUniform()-1.);
	u1 = radius*sin(theta)*cos(phi);
	u2 = radius*sin(theta)*sin(phi);
	u3 = radius*cos(theta);
}
#endif


void PureGaugeUpdater::updateLink(extended_gauge_lattice_t& lattice, int site, int mu, GaugeAction* action, double beta) {
	GaugeGroup staple = action->staple(lattice, site, mu);
#if NUMCOLORS > 2
	//take the plaquette
	GaugeGroup plaquette = lattice[site][mu]*(staple);
	for (unsigned int k = 0; k < numberColors-1; ++k) {
		for (unsigned int l = k+1; l < numberColors; ++l) {
			//Take the su2 subgroup matrix for the cabibbo marinari update
			real_t aeff = (imag(plaquette.at(k,k))-imag(plaquette.at(l,l)))/2.;
			real_t beff = (imag(plaquette.at(k,l))+imag(plaquette.at(l,k)))/2.;
			real_t ceff = (real(plaquette.at(k,l))-real(plaquette.at(l,k)))/2.;
			real_t deff = (real(plaquette.at(k,k))+real(plaquette.at(l,l)))/2.;
			real_t detStaple = aeff*aeff + beff*beff + ceff*ceff + deff*deff;
			real_t b = numberColors/(2.*beta*sqrt(detStaple));
			//Use the standard Kennedy-Pendleton algorithm for su2
			real_t u0, u1, u2, u3;
			u0 = generate_radius(b);
			generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
			//Calculate the su2 update matrix
			matrix2x2_t subupdate;
			subupdate.at(0,0) = std::complex<real_t>(u0, u3);
			subupdate.at(0,1) = std::complex<real_t>(u2, u1);
			subupdate.at(1,0) = std::complex<real_t>(-u2, u1);
			subupdate.at(1,1) = std::complex<real_t>(u0, -u3);
			matrix2x2_t substaple;
			substaple.at(0,0) = std::complex<real_t>(deff, aeff);
			substaple.at(0,1) = std::complex<real_t>(ceff, beff);
			substaple.at(1,0) = std::complex<real_t>(-ceff, beff);
			substaple.at(1,1) = std::complex<real_t>(deff, -aeff);
			//extract the su2 update matrix, compute also one overrelaxation step
			subupdate = htrans(substaple)*htrans((subupdate)*(htrans(substaple)))*htrans(substaple)/pow(detStaple,3./2.);
			//Update the plaquette and the link
			for (unsigned int i = 0; i < numberColors; ++i) {
				std::complex<real_t> tmp1 = subupdate.at(0,0)*lattice[site][mu].at(k,i) + subupdate.at(0,1)*lattice[site][mu].at(l,i);
				std::complex<real_t> tmp2 = subupdate.at(1,0)*lattice[site][mu].at(k,i) + subupdate.at(1,1)*lattice[site][mu].at(l,i);
				lattice[site][mu].at(k,i) = tmp1;
				lattice[site][mu].at(l,i) = tmp2;
				std::complex<real_t> tmp3 = subupdate.at(0,0)*plaquette.at(k,i) + subupdate.at(0,1)*plaquette.at(l,i);
				std::complex<real_t> tmp4 = subupdate.at(1,0)*plaquette.at(k,i) + subupdate.at(1,1)*plaquette.at(l,i);
				plaquette.at(k,i) = tmp3;
				plaquette.at(l,i) = tmp4;
			}
		}
	}
#endif
#if NUMCOLORS == 2
	//compute the parameters for the kennedy algorithm
	real_t detStaple = abs(det(staple));
	real_t b = 1./(beta*sqrt(detStaple));
	real_t u0, u1, u2, u3;
	u0 = generate_radius(b);
	generate_vector(sqrt(1.-u0*u0), u1, u2, u3);
	//update the matrix
	lattice[site][mu].at(0,0) = std::complex<real_t>(u0, u3);
	lattice[site][mu].at(0,1) = std::complex<real_t>(u2, u1);
	lattice[site][mu].at(1,0) = std::complex<real_t>(-u2, u1);
	lattice[site][mu].at(1,1) = std::complex<real_t>(u0, -u3);
	//redefine the matrix U_New = R S^\dag, normalize the determinant
	lattice[site][mu] = (lattice[site][mu])*(htrans(staple))/sqrt(detStaple);
	//perform a single overrelaxation step
	lattice[site][mu] = (htrans(staple)*htrans(lattice[site][mu])*htrans(staple))/detStaple;
#endif
}


} /* namespace Update */
