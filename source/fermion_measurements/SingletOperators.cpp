#include "SingletOperators.h"
#include "io/GlobalOutput.h"
#include "utils/Translate.h"

namespace Update {

SingletOperators::SingletOperators() : StochasticEstimator(), WilsonFlow(), squareDiracOperator(0), diracOperator(0), biConjugateGradient(0) { }

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

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(environment.configurations.get<unsigned int>("generic_inverter_max_steps"));
		biConjugateGradient->setPrecision(environment.configurations.get<real_t>("generic_inverter_precision"));
	}

	std::vector< int* > coordinates;
	{
		int* c = new int[4];
		c[0] = 0;
		c[1] = 0;
		c[2] = 0;
		c[3] = 0;
		coordinates.push_back(c);
	}
	for (int distance = 1; distance < 5; ++distance) {
		int* c = new int[4];
		c[0] = distance;
		c[1] = 0;
		c[2] = 0;
		c[3] = 0;
		coordinates.push_back(c);

		c = new int[4];
		c[0] = 0;
		c[1] = distance;
		c[2] = 0;
		c[3] = 0;
		coordinates.push_back(c);

		c = new int[4];
		c[0] = distance;
		c[1] = distance;
		c[2] = 0;
		c[3] = 0;
		coordinates.push_back(c);

		c = new int[4];
		c[0] = distance;
		c[1] = distance;
		c[2] = distance;
		c[3] = 0;
		coordinates.push_back(c);

		c = new int[4];
		c[0] = distance;
		c[1] = distance;
		c[2] = distance;
		c[3] = distance;
		coordinates.push_back(c);
	}

	long_real_t* A00c = new long_real_t[coordinates.size()];
	long_real_t* A11  = new long_real_t[coordinates.size()];
	long_real_t* A10  = new long_real_t[coordinates.size()];
	long_real_t* A00  = new long_real_t[coordinates.size()];

	for (unsigned int i = 0; i < coordinates.size(); ++i) {
		A00c[i] = 0;
		A11[i]  = 0;
		A10[i]  = 0;
		A00[i]  = 0;
	}

	int inversionSteps = 0;

	extended_dirac_vector_t traceInverse;
	AlgebraUtils::setToZero(traceInverse);
	extended_dirac_vector_t connectedInverse;
	AlgebraUtils::setToZero(connectedInverse);

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise);
		diracOperator->multiply(tmp_square, randomNoise);
		biConjugateGradient->solve(squareDiracOperator, tmp_square, tmp);
		inversionSteps += biConjugateGradient->getLastSteps();
		if (isOutputProcess()) std::cout << "SingletOperators::Inversion " << step << " done in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;
		//diracOperator->multiply(tmp, tmp_square);

		//This part is needed to compute the disconnected contribution
#pragma omp parallel for
		for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					traceInverse[site][mu][c] += conj(randomNoise[site][mu][c])*tmp[site][mu][c];
				}
			}
		}
	}

	//Average
#pragma omp parallel for
	for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				traceInverse[site][mu][c] = traceInverse[site][mu][c]/std::complex<real_t>(max_step,0);
				connectedInverse[site][mu][c] = connectedInverse[site][mu][c]/std::complex<real_t>(max_step,0);
			}
		}
	}

	for (unsigned int index = 0; index < coordinates.size(); ++index) {
		int* coord = coordinates[index];
		extended_dirac_vector_t traceInverseTranslated;
		translate(traceInverse, traceInverseTranslated, coord[0], coord[1], coord[2], coord[3]);
		extended_gauge_lattice_t gaugeTranslated;
		translate(environment.gaugeLinkConfiguration, gaugeTranslated, coord[0], coord[1], coord[2], coord[3]);
	
		long_real_t A11t = 0.;
		long_real_t A10t = 0.;
		long_real_t A00t = 0.;

#pragma omp parallel for reduction(+:A00t,A10t,A11t)
		for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
			std::complex<real_t> disconnectedx = 0;
			std::complex<real_t> disconnectedy = 0;
			//Trace spin color
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					disconnectedx += traceInverse[site][mu][c];
					disconnectedy += traceInverseTranslated[site][mu][c];
				}
			}		
			//Correlation
			A00t += real(disconnectedx)*real(disconnectedy);
					
			A10t += this->measureEnergyAndTopologicalCharge(environment.gaugeLinkConfiguration,site).second*real(disconnectedy);
			A11t += this->measureEnergyAndTopologicalCharge(environment.gaugeLinkConfiguration,site).second*this->measureEnergyAndTopologicalCharge(gaugeTranslated,site).second;
		}
		A00[index] = A00t;
		A10[index] = A10t;
		A11[index] = A11t;
	}
	
	//Now we take the connected part
	extended_dirac_vector_t source, eta;
	extended_dirac_vector_t inverseFull[diracVectorLength*4];
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateSource(source, alpha, c);
			diracOperator->multiply(eta,source);
			biConjugateGradient->solve(squareDiracOperator, eta, inverseFull[c*4 + alpha]);
			//diracOperator->multiply(inverseFull[c*4 + alpha],eta);
			if (isOutputProcess()) std::cout << "SingletOperators::Inversion " << c*4 + alpha << " done in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;
			inversionSteps += biConjugateGradient->getLastSteps();
		}
	}
	
	for (unsigned int index = 0; index < coordinates.size(); ++index) {
		int* coord = coordinates[index];
		//Trace spin color
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			for (int c = 0; c < diracVectorLength; ++c) {
				extended_dirac_vector_t inverseFullTranslated;
				translate(inverseFull[c*4 + alpha], inverseFullTranslated, coord[0], coord[1], coord[2], coord[3]);
				typedef extended_dirac_vector_t::Layout LT;
				
				int site = LT::localIndex[LT::getGlobalCoordinate(0,0,0,0)];
				//Dot
				if (site != -1) {
					for (unsigned int beta = 0; beta < 4; ++beta) {
						for (int d = 0; d < diracVectorLength; ++d) {
							A00c[index] += real(conj(inverseFullTranslated[site][beta][d])*inverseFullTranslated[site][beta][d]);
						}
					}
				}
			}
		}
	}

	

	long_real_t volume = environment.gaugeLinkConfiguration.getLayout().globalVolume;
	for (unsigned int i = 0; i < coordinates.size(); ++i) {
		A00[i] = A00[i]/(volume);
		A10[i] = A10[i]/(volume);
		A11[i] = A11[i]/(volume);
		reduceAllSum(A00c[i]);
		reduceAllSum(A00[i]);
		reduceAllSum(A10[i]);
		reduceAllSum(A11[i]);
	}
	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		
		output->push("pseudoscalar_operator");
		for (unsigned int i = 0; i < coordinates.size(); ++i) {
			int* coord = coordinates[i];
			std::cout << "SingletOperators::Pseudoscalar correlator matrix at x (" << coord[0] << ", " << coord[1] << ", " << coord[2] << ", " << coord[3] << "): {"<< A00[i] + A00c[i] << ", " << A10[i] << ", " << A11[i] << "}" << std::endl;
			std::cout << "SingletOperators::	Connected part of element (0,0): " << A00c[i] << std::endl;
			std::cout << "SingletOperators::	Disconnected part of element (0,0): " << A00[i] << std::endl;

			output->push("pseudoscalar_operator");
			output->write("pseudoscalar_operator", coord[0]);
			output->write("pseudoscalar_operator", coord[1]);
			output->write("pseudoscalar_operator", coord[2]);
			output->write("pseudoscalar_operator", coord[3]);
			output->write("pseudoscalar_operator", A00[i] + A00c[i]);
			output->write("pseudoscalar_operator", A10[i]);
			output->write("pseudoscalar_operator", A11[i]);
			output->pop("pseudoscalar_operator");
		
		}
		output->pop("pseudoscalar_operator");

		output->push("pseudoscalar_operator_connected");
		for (unsigned int i = 0; i < coordinates.size(); ++i) {
			int* coord = coordinates[i];
			output->push("pseudoscalar_operator_connected");
			output->write("pseudoscalar_operator_connected", coord[0]);
			output->write("pseudoscalar_operator_connected", coord[1]);
			output->write("pseudoscalar_operator_connected", coord[2]);
			output->write("pseudoscalar_operator_connected", coord[3]);
			output->write("pseudoscalar_operator_connected", A00c[i]);
			output->write("pseudoscalar_operator_connected", A10[i]);
			output->write("pseudoscalar_operator_connected", A11[i]);
			output->pop("pseudoscalar_operator_connected");
		}
		output->pop("pseudoscalar_operator_connected");
		std::cout << "SingletOperators::Computation done in " << inversionSteps << " inversion steps." << std::endl;
	}

	delete[] A00c;
	delete[] A00;
	delete[] A10;
	delete[] A11;
	for (unsigned int i = 0; i < coordinates.size(); ++i) {
		delete[] coordinates[i];
	}
	
}

}

