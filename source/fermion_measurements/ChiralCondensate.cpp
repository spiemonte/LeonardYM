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

ChiralCondensate::ChiralCondensate() : LatticeSweep(), StochasticEstimator(), squareDiracOperator(0), diracOperator(0), biConjugateGradient(0) { }

ChiralCondensate::ChiralCondensate(const ChiralCondensate& toCopy) : LatticeSweep(toCopy), StochasticEstimator(toCopy), squareDiracOperator(0), diracOperator(0), biConjugateGradient(0) { }

ChiralCondensate::~ChiralCondensate() {
	if (squareDiracOperator) delete squareDiracOperator;
	if (diracOperator) delete diracOperator;
	if (biConjugateGradient) delete biConjugateGradient;
}

void ChiralCondensate::execute(environment_t& environment) {
	unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");

	if (squareDiracOperator == 0) {
		squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	}
	squareDiracOperator->setLattice(environment.getFermionLattice());
	
	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}
	diracOperator->setLattice(environment.getFermionLattice());
	
	long_real_t volume = environment.gaugeLinkConfiguration.getLayout().globalVolume;

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));
	}

	std::vector< long_real_t > chiralCondensateRe;
	std::vector< long_real_t > chiralCondensateIm;
	std::vector< long_real_t > pseudoCondensateRe;
	std::vector< long_real_t > pseudoCondensateIm;
	std::vector< long_real_t > chiralCondensateConnectedRe;
	std::vector< long_real_t > chiralCondensateConnectedIm;

	std::vector< long_real_t > pionNormRe;
	std::vector< long_real_t > pionNormIm;

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise);
		
		biConjugateGradient->solve(squareDiracOperator, randomNoise, tmp_square);
		diracOperator->multiply(tmp, tmp_square);
		
		std::complex<long_real_t> condensate = AlgebraUtils::gamma5dot(randomNoise, tmp);
		chiralCondensateRe.push_back(real(condensate)/volume);
		chiralCondensateIm.push_back(imag(condensate)/volume);
		
		std::complex<long_real_t> pseudoCondensate = AlgebraUtils::dot(randomNoise, tmp);
		pseudoCondensateRe.push_back(real(pseudoCondensate)/volume);
		pseudoCondensateIm.push_back(imag(pseudoCondensate)/volume);

		std::complex<long_real_t> pionNorm = AlgebraUtils::dot(randomNoise, tmp_square);
		pionNormRe.push_back(real(pionNorm)/volume);
		pionNormIm.push_back(imag(pionNorm)/volume);
		
		AlgebraUtils::gamma5(randomNoise);
		biConjugateGradient->solve(squareDiracOperator, randomNoise, tmp_square);
		diracOperator->multiply(randomNoise, tmp_square);

		std::complex<long_real_t> condensateConnected = AlgebraUtils::gamma5dot(randomNoise, tmp);
		chiralCondensateConnectedRe.push_back(real(condensateConnected)/volume);
		chiralCondensateConnectedIm.push_back(imag(condensateConnected)/volume);
	}

	if (isOutputProcess() && environment.measurement) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("condensate");
		output->push("pseudocondensate");
		output->push("condensate_connected");
		output->push("pion_norm");

		std::cout << "ChiralCondensate::Chiral condensate is (re) " << this->mean(chiralCondensateRe) << " +/- " << this->standardDeviation(chiralCondensateRe) << std::endl;
		std::cout << "ChiralCondensate::Chiral condensate is (im) " << this->mean(chiralCondensateIm) << " +/- " << this->standardDeviation(chiralCondensateIm) << std::endl;
		
		std::cout << "ChiralCondensate::Pseudo condensate is (re) " << this->mean(pseudoCondensateRe) << " +/- " << this->standardDeviation(pseudoCondensateRe) << std::endl;
		std::cout << "ChiralCondensate::Pseudo condensate is (im) " << this->mean(pseudoCondensateIm) << " +/- " << this->standardDeviation(pseudoCondensateIm) << std::endl;

		std::cout << "ChiralCondensate::Connected chiral condensate susceptibility is (re) " << this->mean(chiralCondensateConnectedRe) << " +/- " << this->standardDeviation(chiralCondensateConnectedRe) << std::endl;
		std::cout << "ChiralCondensate::Connected chiral condensate susceptibility is (im) " << this->mean(chiralCondensateConnectedIm) << " +/- " << this->standardDeviation(chiralCondensateConnectedIm) << std::endl;

		std::cout << "ChiralCondensate::Pion norm is (re) " << this->mean(pionNormRe) << " +/- " << this->standardDeviation(pionNormRe) << std::endl;
		std::cout << "ChiralCondensate::Pion norm is (im) " << this->mean(pionNormIm) << " +/- " << this->standardDeviation(pionNormIm) << std::endl;

		output->write("condensate", this->mean(chiralCondensateRe));
		output->write("condensate", this->mean(chiralCondensateIm));
		
		output->write("pseudocondensate", this->mean(pseudoCondensateRe));
		output->write("pseudocondensate", this->mean(pseudoCondensateIm));

		output->write("condensate_connected", this->mean(chiralCondensateConnectedRe));
		output->write("condensate_connected", this->mean(chiralCondensateConnectedIm));

		output->write("pion_norm", this->mean(pionNormRe));
		output->write("pion_norm", this->mean(pionNormIm));

		output->pop("condensate");
		output->pop("pseudocondensate");
		output->pop("condensate_connected");
		output->pop("pion_norm");
	}
	
}

} /* namespace Update */
