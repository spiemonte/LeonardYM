#include "SingletOperators.h"
#include "io/GlobalOutput.h"

namespace Update {

SingletOperators::SingletOperators() : StochasticEstimator(), WilsonFlow() { }

void SingletOperators::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t Lt;
	unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");

	if (squareDiracOperator == 0) {
		squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	}
	squareDiracOperator->setLattice(environment.getFermionLattice());
	
	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}
	diracOperator->setLattice(environment.getFermionLattice());

	long_real_t A11 = 0.;
	long_real_t A10 = 0.;
	long_real_t A00 = 0.;

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise);
		biConjugateGradient->solve(squareDiracOperator, randomNoise, tmp_square);
		diracOperator->multiply(tmp, tmp_square);

#pragma omp parallel for reduction(+:A00,A10,A11)
		for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					std::complex<real_t> disconnectedx = conj(randomNoise[site][mu][c])*tmp[site][mu][c];
					std::complex<real_t> disconnectedy = conj(randomNoise[Lt::sup(Lt::sup(Lt::sup(Lt::sup(site,0),1),2),3)][mu][c])*tmp[Lt::sup(Lt::sup(Lt::sup(Lt::sup(site,0),1),2),3)][mu][c];
					std::complex<real_t> connected = conj(randomNoise[Lt::sup(Lt::sup(Lt::sup(Lt::sup(site,0),1),2),3)][mu][c])*tmp[site][mu][c]*randomNoise[site][mu][c]*conj(tmp[Lt::sup(Lt::sup(Lt::sup(Lt::sup(site,0),1),2),3)][mu][c]);
					A00 += real(disconnectedx*disconnectedy)+real(connected);
					
					A10 += this->measureEnergyAndTopologicalCharge(environment.gaugeLinkConfiguration,site).second*real(disconnectedy);
					A11 += this->measureEnergyAndTopologicalCharge(environment.gaugeLinkConfiguration,site).second*this->measureEnergyAndTopologicalCharge(environment.gaugeLinkConfiguration,Lt::sup(Lt::sup(Lt::sup(Lt::sup(site,0),1),2),3)).second;
				}
			}
		}
	}


	long_real_t volume = environment.gaugeLinkConfiguration.getLayout().globalVolume;
	A00 = A00/(volume*max_step);
	A10 = A10/(volume*max_step);
	A11 = A11/(volume*max_step);

	std::cout << "Vediamo di che morte dobbiamo morire: " << A00 << " " << A10 << " " << A11 << std::endl;
}

}

