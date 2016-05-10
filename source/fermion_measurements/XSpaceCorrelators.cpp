#include "XSpaceCorrelators.h"
#include "io/GlobalOutput.h"
#include "utils/Translate.h"

namespace Update {

XSpaceCorrelators::XSpaceCorrelators() : StochasticEstimator(), WilsonFlow(), squareDiracOperator(0), diracOperator(0), biConjugateGradient(0), gammaOperators() { }

void XSpaceCorrelators::execute(environment_t& environment) {
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
	for (int distance = 1; distance < 4; ++distance) {
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
		c[0] = 0;
		c[1] = 0;
		c[2] = 0;
		c[3] = distance;
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

	int inversionSteps = 0;

	extended_dirac_vector_t traceInverse;
	AlgebraUtils::setToZero(traceInverse);
	extended_dirac_vector_t* randomNoise = new extended_dirac_vector_t[max_step];
	extended_dirac_vector_t* inverseRandomNoise = new extended_dirac_vector_t[max_step];

	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);

	for (unsigned int step = 0; step < max_step; ++step) {
		this->generateRandomNoise(randomNoise[step]);
		tmp = randomNoise[step];
		AlgebraUtils::gamma5(tmp);
		diracOperator->multiply(source, tmp);
		biConjugateGradient->solve(squareDiracOperator, source, inverseRandomNoise[step]);
		
		inversionSteps += biConjugateGradient->getLastSteps();
		if (isOutputProcess()) std::cout << "XSpaceCorrelators::Inversion " << step << " done in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;

		//This part is needed to compute the disconnected contribution
#pragma omp parallel for
		for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (int c = 0; c < diracVectorLength; ++c) {
					traceInverse[site][mu][c] += conj(randomNoise[step][site][mu][c])*inverseRandomNoise[step][site][mu][c];
				}
			}
		}
	}

#pragma omp parallel for
	for (int site = 0; site < environment.getFermionLattice().completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (int c = 0; c < diracVectorLength; ++c) {
				traceInverse[site][mu][c] = traceInverse[site][mu][c]/static_cast<real_t>(max_step);
			}
		}
	}

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	if (isOutputProcess()) std::cout << "XSpaceCorrelators::Inversions done in: " << (elapsed) << " s."<< std::endl;
	if (isOutputProcess()) std::cout << "XSpaceCorrelators::Computation done in " << inversionSteps << " inversion steps." << std::endl;

	clock_gettime(CLOCK_REALTIME, &start);

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		for (int gamma = 0; gamma < 16; ++gamma) {
			std::string name = "meson_correlator_";
			name += toString(gamma);

			output->push(name);

			name = "meson_correlator_disconnected_";
			name += toString(gamma);

			output->push(name);
		}
	}
	
	for (unsigned int index = 0; index < coordinates.size(); ++index) {
		int* coord = coordinates[index];

		std::complex<long_real_t> result_connected[16][max_step];
		std::complex<long_real_t> result_disconnected[16];

		for (unsigned int step = 0; step < max_step; ++step) {
			extended_dirac_vector_t randomNoiseTranslated, inverseRandomNoiseTranslated, propagator, propagatorInverse;
			translate(randomNoise[step], randomNoiseTranslated, coord[0], coord[1], coord[2], coord[3]);
			translate(inverseRandomNoise[step], inverseRandomNoiseTranslated, coord[0], coord[1], coord[2], coord[3]);

/*#pragma omp parallel for
			for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						propagator[site][mu][c] = conj(randomNoiseTranslated[site][mu][c])*inverseRandomNoise[step][site][mu][c];
						propagatorInverse[site][mu][c] = conj(randomNoise[step][site][mu][c])*inverseRandomNoiseTranslated[site][mu][c];
					}
				}
			}*/


			for (int gamma = 0; gamma < 16; ++gamma) {
				extended_dirac_vector_t /*gammaPropagator, gammaPropagatorInverse,*/ hRandomNoise, gammaRandomNoise, gammaRandomNoiseTranslated;
				AlgebraUtils::conjugate(hRandomNoise, randomNoiseTranslated);
				gammaOperators.multiply(gammaRandomNoiseTranslated, hRandomNoise, gamma);
				
				AlgebraUtils::conjugate(hRandomNoise, randomNoise[step]);
				gammaOperators.multiply(gammaRandomNoise, hRandomNoise, gamma);

				matrix_t gr   = gammas.gammaChromaMatrices(gamma);
				matrix_t grh = gammas.gammaChromaMatrices(8)*htrans(gr)*gammas.gammaChromaMatrices(8);

				long_real_t tmp_re = 0, tmp_im = 0;
#pragma omp parallel for reduction(+:tmp_re, tmp_im)
				for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
					std::complex<real_t> tmp(0.);
					for (unsigned int mu = 0; mu < 4; ++mu) {
						for (unsigned int nu = 0; nu < 4; ++nu) {
							if (gr.at(mu,nu) != std::complex<real_t>(0.,0.)) {
								for (unsigned int alpha = 0; alpha < 4; ++alpha) {
									for (unsigned int beta = 0; beta < 4; ++beta) {
										if (grh.at(alpha,beta) != std::complex<real_t>(0.,0.)) {
											for (int c1 = 0; c1 < diracVectorLength; ++c1) {
												for (int c2 = 0; c2 < diracVectorLength; ++c2) {
									//tmp +=  gammaRandomNoiseTranslated[site][mu][c1]*inverseRandomNoise[step][site][nu][c2]*gammaRandomNoise[site][nu][c2]*inverseRandomNoiseTranslated[site][mu][c1];
													tmp += gr.at(mu,nu)*conj(randomNoiseTranslated[site][alpha][c2])*inverseRandomNoise[step][site][nu][c1]*grh.at(alpha,beta)*conj(randomNoise[step][site][mu][c1])*inverseRandomNoiseTranslated[site][beta][c2];
												}
											}
										}
									}
								}
							}
						}
					}
					tmp_re += tmp.real();
					tmp_im += tmp.imag();
				}
				reduceAllSum(tmp_re);
				reduceAllSum(tmp_im);

				result_connected[gamma][step] = std::complex<long_real_t>(tmp_re,tmp_im);
			}
		}
		
		extended_dirac_vector_t traceInverseTranslated;
		translate(traceInverse, traceInverseTranslated, coord[0], coord[1], coord[2], coord[3]);
		for (int gamma = 0; gamma < 16; ++gamma) {
			matrix_t gr   = gammas.gammaChromaMatrices(gamma);
			matrix_t grh = gammas.gammaChromaMatrices(8)*htrans(gr)*gammas.gammaChromaMatrices(8);

			extended_dirac_vector_t gammaTraceInverse, gammaTraceInverseTranslated;
			gammaOperators.multiply(gammaTraceInverse, traceInverse, gamma);
			gammaOperators.multiply(gammaTraceInverseTranslated, traceInverseTranslated, gamma);

			long_real_t tmp_re = 0, tmp_im = 0;
#pragma omp parallel for reduction(+:tmp_re, tmp_im)
			for (int site = 0; site < environment.getFermionLattice().localsize; ++site) {
				std::complex<real_t> disconnectedx = 0;
				std::complex<real_t> disconnectedy = 0;
				//Trace spin color
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (int c = 0; c < diracVectorLength; ++c) {
						disconnectedx += gammaTraceInverse[site][mu][c];
						disconnectedy += gammaTraceInverseTranslated[site][mu][c];
					}
				}		
				//Correlation
				tmp_re += real(disconnectedx*disconnectedy);
				tmp_im += imag(disconnectedx*disconnectedy);
			}
			reduceAllSum(tmp_re);
			reduceAllSum(tmp_im);

			result_disconnected[gamma] = std::complex<long_real_t>(tmp_re,tmp_im);
		}


		if (environment.measurement && isOutputProcess()) {
			//Correct normalization
			long_real_t factor = 4*diracOperator->getKappa()*diracOperator->getKappa();
			GlobalOutput* output = GlobalOutput::getInstance();

			for (int gamma = 0; gamma < 16; ++gamma) {
				std::string name = "meson_correlator_";
				name += toString(gamma);

				output->push(name);

				output->push(name);
				output->write(name, coord[0]);
				output->write(name, coord[1]);
				output->write(name, coord[2]);
				output->write(name, coord[3]);
				output->pop(name);

				for (unsigned int step = 0; step < max_step; ++step) {
					output->push(name);
					output->write(name, factor*result_connected[gamma][step].real());
					output->write(name, factor*result_connected[gamma][step].imag());
					output->pop(name);
				}
				output->pop(name);

				name = "meson_correlator_disconnected_";
				name += toString(gamma);

				output->push(name);

				output->push(name);
				output->write(name, coord[0]);
				output->write(name, coord[1]);
				output->write(name, coord[2]);
				output->write(name, coord[3]);
				output->pop(name);

				output->push(name);
				output->write(name, factor*result_disconnected[gamma].real());
				output->write(name, factor*result_disconnected[gamma].imag());
				output->pop(name);
				
				output->pop(name);
			}
		
		}
		
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		
		for (int gamma = 0; gamma < 16; ++gamma) {
			std::string name = "meson_correlator_";
			name += toString(gamma);

			output->pop(name);

			name = "meson_correlator_disconnected_";
			name += toString(gamma);

			output->pop(name);
		}
	}

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	if (isOutputProcess()) std::cout << "XSpaceCorrelators::Contractions done in: " << (elapsed) << " s."<< std::endl;

	//Now we take the connected part
	extended_dirac_vector_t source, eta;
	extended_dirac_vector_t inverseFull[diracVectorLength*4];
	
	for (unsigned int alpha = 0; alpha < 4; ++alpha) {
		for (int c = 0; c < diracVectorLength; ++c) {
			this->generateSource(source, alpha, c);
			tmp = source;
			AlgebraUtils::gamma5(tmp);
			diracOperator->multiply(source, tmp);
			biConjugateGradient->solve(squareDiracOperator, source, inverseFull[c*4 + alpha]);
			

			//diracOperator->multiply(eta,source);
			//biConjugateGradient->solve(squareDiracOperator, eta, inverseFull[c*4 + alpha]);
			//AlgebraUtils::gamma5(inverseFull[c*4 + alpha]);
			//diracOperator->multiply(inverseFull[c*4 + alpha],eta);
			if (isOutputProcess()) std::cout << "XSpaceCorrelators::Inversion " << c*4 + alpha << " done in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;
			inversionSteps += biConjugateGradient->getLastSteps();
		}
	}
	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		for (int gamma = 0; gamma < 16; ++gamma) {
			std::string name = "meson_correlator_exact_";
			name += toString(gamma);

			output->push(name);
		}
	}

	for (unsigned int index = 0; index < coordinates.size(); ++index) {
		int* coord = coordinates[index];

		std::complex<long_real_t> result_connected[16];
		
		for (int gamma_index = 0; gamma_index < 16; ++gamma_index) {
			result_connected[gamma_index] = 0.;
			//Trace spin color
			for (unsigned int rho = 0; rho < 4; ++rho) {
				for (int c = 0; c < diracVectorLength; ++c) {
					//extended_dirac_vector_t inverseFullTranslated;
					//translate(inverseFull[c*4 + alpha], inverseFullTranslated, coord[0], coord[1], coord[2], coord[3]);
					typedef extended_dirac_vector_t::Layout LT;
				
					int site = LT::localIndex[LT::getGlobalCoordinate(coord[0], coord[1], coord[2], coord[3])];
					//Dot
					if (site != -1) {
						for (unsigned int alpha = 0; alpha < 4; ++alpha) {
							for (unsigned int mu = 0; mu < 4; ++mu) {
								for (unsigned int nu = 0; nu < 4; ++nu) {
									for (int d = 0; d < diracVectorLength; ++d) {
										real_t sign = 1.;
										if (alpha < 2) sign = -sign;
										if (mu < 2) sign = -sign;
										
										result_connected[gamma_index] += sign*gammas.gammaChromaMatrices(gamma_index).at(rho, alpha)*conj(inverseFull[c*4+alpha][site][mu][d])*gammas.gammaChromaMatrices(gamma_index).at(mu, nu)*inverseFull[d*4+nu][site][rho][c];
									}
								}
							}
						}
					}
				}
			}
			reduceAllSum(result_connected[gamma_index].real());
			reduceAllSum(result_connected[gamma_index].imag());
		}

		if (environment.measurement && isOutputProcess()) {
			GlobalOutput* output = GlobalOutput::getInstance();
			long_real_t factor = 4*diracOperator->getKappa()*diracOperator->getKappa();

			for (int gamma = 0; gamma < 16; ++gamma) {
				std::string name = "meson_correlator_exact_";
				name += toString(gamma);

				output->push(name);

				output->push(name);
				output->write(name, coord[0]);
				output->write(name, coord[1]);
				output->write(name, coord[2]);
				output->write(name, coord[3]);
				output->pop(name);

				output->push(name);
				output->write(name, factor*result_connected[gamma].real());
				output->write(name, factor*result_connected[gamma].imag());
				output->pop(name);
				
				output->pop(name);
			}
		}

		
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		for (int gamma = 0; gamma < 16; ++gamma) {
			std::string name = "meson_correlator_exact_";
			name += toString(gamma);

			output->pop(name);
		}
	}


	
	for (unsigned int i = 0; i < coordinates.size(); ++i) {
		delete[] coordinates[i];
	}
	delete[] randomNoise;
	delete[] inverseRandomNoise;
	
}

}

