#include "LandauGluonPropagator.h"
#include "ToString.h"
#include "io/GlobalOutput.h"
#include "utils/LieGenerators.h"

namespace Update {

LandauGluonPropagator::LandauGluonPropagator() : LandauGaugeFixing() { }

LandauGluonPropagator::LandauGluonPropagator(const LandauGluonPropagator& toCopy) : LandauGaugeFixing(toCopy) { }

LandauGluonPropagator::~LandauGluonPropagator() { }

void LandauGluonPropagator::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;//TODO: only vector operations?
	typedef extended_gauge_lattice_t LT;

	std::vector<real_t> maximal_momentum = environment.configurations.get< std::vector<real_t> >("LandauGluonPropagator::maximal_momentum");
	real_t cut =  environment.configurations.get< real_t >("LandauGluonPropagator::anisotropy_cut");

	std::vector< std::vector<real_t> > momenta;
	for (int i = 0; i < maximal_momentum[0]; ++i) {
		for (int j = 0; j < maximal_momentum[1]; ++j) {
			for (int k = 0; k < maximal_momentum[2]; ++k) {
				for (int l = 0; l < maximal_momentum[3]; ++l) {
					std::vector<real_t> momentum(4);
					momentum[0] = i*2.*PI/Layout::glob_x;
					momentum[1] = j*2.*PI/Layout::glob_y;
					momentum[2] = k*2.*PI/Layout::glob_z;
					momentum[3] = l*2.*PI/Layout::glob_t;

					real_t mean = (momentum[0] + momentum[1] + momentum[2] + momentum[3])/4.;
					real_t deviation = fabs(momentum[0] - mean) + fabs(momentum[1] - mean) + fabs(momentum[2] - mean) + fabs(momentum[3] - mean);
					if (deviation < cut) {
						momenta.push_back(momentum);
					}
				}
			}
		}
	}

	std::vector< std::vector<real_t> > delta(4);
	for (int mu = 0; mu < 4; ++mu) {
		std::vector<real_t> delta_mu(4);
		for (int nu = 0; nu < 4; ++nu) {
			if (mu == nu) {
				delta_mu[nu] = 1./2.;
			} 
			else {
				delta_mu[nu] = 0.;
			}
		}
		delta[mu] = delta_mu;
	}

	real_t convergence = this->deviation(environment.gaugeLinkConfiguration);

	if (isOutputProcess()) std::cout << "LandauGluonPropagator::Deviation from the Landau gauge: " << convergence << std::endl;
	if (isOutputProcess()) std::cout << "LandauGluonPropagator::Analizing " << momenta.size() << " momenta configurations ..." << std::endl;
	

	extended_gauge_lattice_t Afield;
	this->getLieAlgebraField(Afield, environment.gaugeLinkConfiguration);


	LieGenerator<GaugeGroup> lieGenerators;

	std::complex<real_t> Azero[4][lieGenerators.numberGenerators()];

	std::vector< std::vector<real_t> > data;

	for (unsigned int i = 0; i < momenta.size(); ++i) {
		std::vector<real_t> momentum = momenta[i];

		std::complex<real_t> resultAk[4][lieGenerators.numberGenerators()], resultAmk[4][lieGenerators.numberGenerators()];
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
				resultAk[mu][c] = 0.;
				resultAmk[mu][c] = 0.;
			}
		}

#ifdef MULTITHREADING
		std::complex<real_t> Ak[omp_get_max_threads()][4][lieGenerators.numberGenerators()], Amk[omp_get_max_threads()][4][lieGenerators.numberGenerators()];
		for (int k = 0; k < omp_get_max_threads(); ++k) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
					Ak[k][mu][c] = 0.;
					Amk[k][mu][c] = 0.;
				}
			}
		}
#endif

#pragma omp parallel for
		for (int site = 0; site < Afield.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				real_t phase = 
					   (Layout::globalIndexX(site) + delta[mu][0])*momentum[0]
					+ (Layout::globalIndexY(site) + delta[mu][1])*momentum[1] 
					+ (Layout::globalIndexZ(site) + delta[mu][2])*momentum[2] 
					+ (Layout::globalIndexT(site) + delta[mu][3])*momentum[3] ;
				for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
#ifndef MULTITHREADING
					resultAk[mu][c] += std::complex<real_t>(cos(phase),sin(phase))*trace(Afield[site][mu]*lieGenerators.get(c));
					resultAmk[mu][c] += std::complex<real_t>(cos(phase),-sin(phase))*trace(Afield[site][mu]*lieGenerators.get(c));
#endif
#ifdef MULTITHREADING
					Ak[omp_get_thread_num()][mu][c] += std::complex<real_t>(cos(phase),sin(phase))*trace(Afield[site][mu]*lieGenerators.get(c));
					Amk[omp_get_thread_num()][mu][c] += std::complex<real_t>(cos(phase),-sin(phase))*trace(Afield[site][mu]*lieGenerators.get(c));
#endif
				}
			}
		}

		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
#ifdef MULTITHREADING
				for (int k = 0; k < omp_get_max_threads(); ++k) {
					resultAk[mu][c] += Ak[k][mu][c];
					resultAmk[mu][c] += Amk[k][mu][c];
				}
#endif
				reduceAllSum(resultAk[mu][c]);
				reduceAllSum(resultAmk[mu][c]);
			}
		}

		//Physical momentum
		std::vector<real_t> p(4);
		p[0] = 2.*sin(momentum[0]/2.);
		p[1] = 2.*sin(momentum[1]/2.);
		p[2] = 2.*sin(momentum[2]/2.);
		p[3] = 2.*sin(momentum[3]/2.);

		real_t modulus = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);

		std::vector<real_t> result;
		result.push_back(modulus);

		real_t propagator = 0.;
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
				propagator += real(resultAk[mu][c]*resultAmk[mu][c]);
			}
		}

		result.push_back(propagator);
	
		//The normalization at zero momentum is 4
		if (modulus < 0.0000000000000001) {
			result[1] = result[1]/4.;
		}
		//Otherwise 3
		else {
			result[1] = result[1]/3.;
		}

		if (isOutputProcess()) std::cout << "Momenta p = {" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "} |p| = " << result[0] << std::endl;
		if (isOutputProcess()) std::cout << "                     G2(p,-p) = " << result[1] << std::endl;

		if (i == 0) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < lieGenerators.numberGenerators(); ++c) {
					Azero[mu][c] = resultAk[mu][c];
				}
			}
			result.push_back(0);
		}
		else {
			std::complex<real_t> three_point_function = 0.;
			for (unsigned int c1 = 0; c1 < lieGenerators.numberGenerators(); ++c1) {
				for (unsigned int c2 = 0; c2 < lieGenerators.numberGenerators(); ++c2) {
					for (unsigned int c3 = 0; c3 < lieGenerators.numberGenerators(); ++c3) {
						for (unsigned int mu1 = 0; mu1 < 4; ++mu1) {
							for (unsigned int mu2 = 0; mu2 < 4; ++mu2) {
								three_point_function += std::complex<real_t>(0.,-2. ) * trace((lieGenerators.get(c1)*lieGenerators.get(c2) - lieGenerators.get(c2)*lieGenerators.get(c1))*lieGenerators.get(c3)) * resultAk[mu1][c1] * Azero[mu2][c2] * resultAmk[mu1][c3] * p[mu2] / (6.*modulus*modulus);
								/*std::complex<real_t> che_cazzo = ( std::complex<real_t>(0.,-2.) * trace((lieGenerators.get(c1)*lieGenerators.get(c2) - lieGenerators.get(c2)*lieGenerators.get(c1))*lieGenerators.get(c3)) * resultAk[mu1][c1] * Azero[mu2][c2] * resultAmk[mu1][c3] * p[mu2] / (6.*modulus*modulus) );
								if (std::abs(che_cazzo) > 0.000000000001) {
									std::cout << "Giusto per: " << mu2 << " " << mu1 << " " << c1 << " " << c2 << " " << c3 << "  -  " << p[mu2] << "   " << trace((lieGenerators.get(c1)*lieGenerators.get(c2) - lieGenerators.get(c2)*lieGenerators.get(c1))*lieGenerators.get(c3)) << "   " << resultAk[mu1][c1] << "   " << Azero[mu2][c2]  << "   " <<  resultAmk[mu1][c3] << std::endl;
									std::cout << " E poi: " << che_cazzo << " " << three_point_function << std::endl << std::endl;
								}*/
							}
						}
					}
				}
			}

			/*std::complex<real_t> test1 = 0., test2 = 0.;
			for (unsigned int c1 = 0; c1 < lieGenerators.numberGenerators(); ++c1) {
				for (unsigned int c2 = 0; c2 < lieGenerators.numberGenerators(); ++c2) {
					for (unsigned int c3 = 0; c3 < lieGenerators.numberGenerators(); ++c3) {
						test1 += 2. * trace((lieGenerators.get(c1)*lieGenerators.get(c2) - lieGenerators.get(c2)*lieGenerators.get(c1))*lieGenerators.get(c3)) * resultAk[0][c1] * Azero[1][c2] * resultAmk[0][c3] * p[1] / (6.*modulus*modulus);
						test2 += 2. * trace((lieGenerators.get(c1)*lieGenerators.get(c2) - lieGenerators.get(c2)*lieGenerators.get(c1))*lieGenerators.get(c3)) * resultAk[1][c1] * Azero[1][c2] * resultAmk[1][c3] * p[1] / (6.*modulus*modulus);
					}
				}
			}

			std::cout << "                 Test: " << test1/test2 << std::endl;*/
			
			result.push_back(imag(three_point_function));
			if (isOutputProcess()) std::cout << "                     G3(p,0,-p) = " << three_point_function << std::endl;
		}

		data.push_back(result);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();
		output->push("gluon_propagator");

		for (unsigned int i = 0; i < data.size(); ++i) {
			bool already_considered = false;
			for (unsigned int  j = 0; j < i; ++j) {
				if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
					already_considered = true;
					break;
				}
			}
			if (!already_considered) {
				std::cout << "{";
				int nn = 1;
				real_t averaged_result = data[i][1];
				
				std::vector<real_t> p(4);
				p[0] = sin(momenta[i][0]/2.);
				p[1] = sin(momenta[i][1]/2.);
				p[2] = sin(momenta[i][2]/2.);
				p[3] = sin(momenta[i][3]/2.);
				real_t normalization = 4.*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3] + (1./3.)*(p[0]*p[0]*p[0]*p[0] + p[1]*p[1]*p[0]*p[0] + p[2]*p[2]*p[2]*p[2] + p[3]*p[3]*p[3]*p[3]));
				std::cout << "{" << momenta[i][0] << ", " << momenta[i][1] << ", " << momenta[i][2] << ", " << momenta[i][3] << ", " << normalization << "},";

				for (unsigned int  j = i + 1; j < data.size(); ++j) {
					if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
						p[0] = sin(momenta[j][0]/2.);
						p[1] = sin(momenta[j][1]/2.);
						p[2] = sin(momenta[j][2]/2.);
						p[3] = sin(momenta[j][3]/2.);
						normalization = 4.*(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3] + (1./3.)*(p[0]*p[0]*p[0]*p[0] + p[1]*p[1]*p[0]*p[0] + p[2]*p[2]*p[2]*p[2] + p[3]*p[3]*p[3]*p[3]));
						std::cout << "{" << momenta[j][0] << ", " << momenta[j][1] << ", " << momenta[j][2] << ", " << momenta[j][3] << ", " << normalization << "},";
						++nn;
						averaged_result += data[j][1];
					}
				}
				std::cout << "}," << std::endl;
				output->push("gluon_propagator");
				output->write("gluon_propagator", data[i][0]);
				output->write("gluon_propagator", averaged_result/nn);
				output->pop("gluon_propagator");
			}
		}
		
		output->pop("gluon_propagator");



		output->push("three_gluon_vertex");

		for (unsigned int i = 0; i < data.size(); ++i) {
			bool already_considered = false;
			for (unsigned int  j = 0; j < i; ++j) {
				if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
					already_considered = true;
					break;
				}
			}
			if (!already_considered) {
				int nn = 1;
				real_t averaged_result = data[i][2];
				for (unsigned int  j = i + 1; j < data.size(); ++j) {
					if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
						++nn;
						averaged_result += data[j][2];
					}
				}
				output->push("three_gluon_vertex");
				output->write("three_gluon_vertex", data[i][0]);
				output->write("three_gluon_vertex", averaged_result/nn);
				output->pop("three_gluon_vertex");
			}
		}

		output->pop("three_gluon_vertex");
	}

}

void LandauGluonPropagator::registerParameters(po::options_description& desc) {
	desc.add_options()
		("LandauGluonPropagator::anisotropy_cut", po::value<real_t>()->default_value(0.25), "Cut for the off-diagonal momenta")
		("LandauGluonPropagator::maximal_momentum", po::value<std::string>()->default_value("{2,2,2,2}"), "Momentum for the measure of the vertex function (syntax: {px,py,pz,pt})")
		;
}

}

