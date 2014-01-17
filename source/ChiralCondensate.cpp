/*
 * ChiralCondensate.cpp
 *
 *  Created on: Jul 23, 2012
 *      Author: spiem_01
 */

#include "ChiralCondensate.h"
#include "AlgebraUtils.h"
#include "BiConjugateGradient.h"
#include "GlobalOutput.h"

namespace Update {

ChiralCondensate::ChiralCondensate() : LatticeSweep(), StochasticEstimator(), diracOperator(0), biConjugateGradient(0) { }

ChiralCondensate::ChiralCondensate(const ChiralCondensate& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), diracOperator(0), biConjugateGradient(0) { }

ChiralCondensate::~ChiralCondensate() {
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
}

void ChiralCondensate::execute(environment_t& environment) {
	unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");

	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	diracOperator->setLattice(environment.getFermionLattice());
	diracOperator->setGamma5(false);
	
	long_real_t volume = environment.gaugeLinkConfiguration.getLayout().globalVolume;

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));
	}

	std::vector< long_real_t > chiralCondensateRe;
	std::vector< long_real_t > chiralCondensateIm;
	std::vector< long_real_t > chiralCondensateConnectedRe;
	std::vector< long_real_t > chiralCondensateConnectedIm;

	std::vector< long_real_t > pionNormRe;
	std::vector< long_real_t > pionNormIm;

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise);
		
		biConjugateGradient->solve(diracOperator, randomNoise, tmp);
		biConjugateGradient->solve(diracOperator, tmp, tmp_square);
		
		std::complex<long_real_t> condensate = AlgebraUtils::dot(randomNoise, tmp);
		chiralCondensateRe.push_back(real(condensate)/volume);
		chiralCondensateIm.push_back(imag(condensate)/volume);

		std::complex<long_real_t> condensateConnected = AlgebraUtils::dot(randomNoise, tmp_square);
		chiralCondensateConnectedRe.push_back(real(condensateConnected)/volume);
		chiralCondensateConnectedIm.push_back(imag(condensateConnected)/volume);

		std::complex<long_real_t> pionNorm = AlgebraUtils::dot(tmp, tmp);
		pionNormRe.push_back(real(pionNorm)/volume);
		pionNormIm.push_back(imag(pionNorm)/volume);

	}

	if (isOutputProcess() && environment.measurement) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("condensate");
		output->push("condensate_connected");
		output->push("pion_norm");

		std::cout << "Chiral condensate is (re) " << this->mean(chiralCondensateRe) << " +/- " << this->standardDeviation(chiralCondensateRe) << std::endl;
		std::cout << "Chiral condensate is (im) " << this->mean(chiralCondensateIm) << " +/- " << this->standardDeviation(chiralCondensateIm) << std::endl;

		std::cout << "Connected chiral condensate susceptibility is (re) " << this->mean(chiralCondensateConnectedRe) << " +/- " << this->standardDeviation(chiralCondensateConnectedRe) << std::endl;
		std::cout << "Connected chiral condensate susceptibility is (im) " << this->mean(chiralCondensateConnectedIm) << " +/- " << this->standardDeviation(chiralCondensateConnectedIm) << std::endl;

		std::cout << "Pion norm is (re) " << this->mean(pionNormRe) << " +/- " << this->standardDeviation(pionNormRe) << std::endl;
		std::cout << "Pion norm is (im) " << this->mean(pionNormIm) << " +/- " << this->standardDeviation(pionNormIm) << std::endl;

		output->write("condensate", this->mean(chiralCondensateRe));
		output->write("condensate", this->mean(chiralCondensateIm));

		output->write("condensate_connected", this->mean(chiralCondensateConnectedRe));
		output->write("condensate_connected", this->mean(chiralCondensateConnectedIm));

		output->write("pion_norm", this->mean(pionNormRe));
		output->write("pion_norm", this->mean(pionNormIm));

		output->pop("condensate");
		output->pop("condensate_connected");
		output->pop("pion_norm");
	}
	
}

} /* namespace Update */
