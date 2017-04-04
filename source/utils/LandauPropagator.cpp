#include "LandauPropagator.h"
#include "ToString.h"
#include "io/GlobalOutput.h"
#include "utils/LieGenerators.h"

namespace Update {

LandauPropagator::LandauPropagator() : LandauGaugeFixing() { }

LandauPropagator::LandauPropagator(const LandauPropagator& toCopy) : LandauGaugeFixing(toCopy) { }

LandauPropagator::~LandauPropagator() { }

void LandauPropagator::ghostMatrix(extended_adjoint_color_vector_t& output, const extended_adjoint_color_vector_t& input, const extended_adjoint_lattice_t& A, const extended_adjoint_lattice_t& B, const extended_adjoint_lattice_t& C) const {
	typedef extended_adjoint_color_vector_t LT;

#pragma omp parallel for
	for (int site = 0; site < input.localsize; ++site) {
		set_to_zero(output[site]);
		for (unsigned int mu = 0; mu < 4; ++mu) {
			output[site] += A[site][mu]*input[site] - B[site][mu]*input[LT::sup(site,mu)] - C[site][mu]*input[LT::sdn(site,mu)];
		}
	}

	output.updateHalo();
}

bool LandauPropagator::ghostPropagator(std::complex<real_t>& result, const extended_gauge_lattice_t& lattice, const std::vector<real_t>& momentum, int c, real_t epsilon, unsigned int max_steps) const {
	extended_adjoint_lattice_t A, B, C;
	typedef extended_adjoint_color_vector_t LT;
	typedef extended_adjoint_color_vector_t::Layout Layout;

	LieGenerator<GaugeGroup> lieGenerators;
	
#pragma omp parallel for
	for (int site = 0; site < A.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
				for (unsigned int b = 0; b < lieGenerators.numberGenerators(); ++b) {
					A[site][mu].at(a,b) = real(trace( (lieGenerators.get(a)*lieGenerators.get(b) + lieGenerators.get(b)*lieGenerators.get(a)) *(lattice[site][mu] + lattice[LT::sdn(site,mu)][mu]) ));
					B[site][mu].at(a,b) = 2.*real(trace( lieGenerators.get(b)*lieGenerators.get(a)*lattice[site][mu] ));
					C[site][mu].at(a,b) = 2.*real(trace( lieGenerators.get(a)*lieGenerators.get(b)*lattice[LT::sdn(site,mu)][mu] ));
				}
			}
		}
	}
	A.updateHalo();
	B.updateHalo();
	C.updateHalo();

	extended_adjoint_color_vector_t ghost;
	extended_adjoint_color_vector_t source, tmp, r, p;
#pragma omp parallel for
	for (int site = 0; site < source.localsize; ++site) {
		real_t phase = 	   Layout::globalIndexX(site)*momentum[0]
					+ Layout::globalIndexY(site)*momentum[1] 
					+ Layout::globalIndexZ(site)*momentum[2] 
					+ Layout::globalIndexT(site)*momentum[3] ;
		set_to_zero(source[site]);
		source[site][c] = std::complex<real_t>(cos(phase),sin(phase));
	}

	ghost = source;
	
	this->ghostMatrix(tmp,ghost,A,B,C);

#pragma omp parallel for
	for (int site = 0; site < source.completesize; ++site) {
		r[site] = source[site] - tmp[site];
		p[site] = r[site];
	}

	long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
	for (int site = 0; site < source.localsize; ++site) {
		for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
			norm += real(conj(r[site][a])*r[site][a]);
		}
	}
	reduceAllSum(norm);
	
	long_real_t norm_next = norm;

	for (unsigned int step = 0; step < max_steps; ++step) {
		this->ghostMatrix(tmp,p,A,B,C);
		norm = norm_next;
		long_real_t ddot_re = 0., ddot_im = 0.;
#pragma omp parallel for reduction(+:ddot_re, ddot_im)
		for (int site = 0; site < source.localsize; ++site) {
			for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
				ddot_re += real(conj(p[site][a])*tmp[site][a]);
				ddot_im += imag(conj(p[site][a])*tmp[site][a]);
			}
		}
		reduceAllSum(ddot_re);
		reduceAllSum(ddot_im);
		std::complex<real_t> alpha = std::complex<real_t>(norm,0.)/std::complex<real_t>(ddot_re,ddot_im);


#pragma omp parallel for
		for (int site = 0; site < source.completesize; ++site) {
			ghost[site] = ghost[site] + alpha * p[site];
			r[site] = r[site] - alpha * tmp[site];
		}

		norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
		for (int site = 0; site < r.localsize; ++site) {
			for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
				norm_next += real(conj(r[site][a])*r[site][a]);
			}
		}
		reduceAllSum(norm_next);
		
		if (norm_next < epsilon) {
			long_real_t result_re = 0., result_im = 0.;			
#pragma omp parallel for reduction(+:result_re,result_im)
			for (int site = 0; site < r.localsize; ++site) {
				real_t phase = 	   Layout::globalIndexX(site)*momentum[0]
					+ Layout::globalIndexY(site)*momentum[1] 
					+ Layout::globalIndexZ(site)*momentum[2] 
					+ Layout::globalIndexT(site)*momentum[3] ;
				
				result_re += real(std::complex<real_t>(cos(phase),-sin(phase))*ghost[site][c]);
				result_im += imag(std::complex<real_t>(cos(phase),-sin(phase))*ghost[site][c]);
			}
			reduceAllSum(result_re);
			reduceAllSum(result_im);
			result = std::complex<real_t>(result_re,result_im);


			this->ghostMatrix(tmp,ghost,A,B,C);
			norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
			for (int site = 0; site < r.localsize; ++site) {
				for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
					norm_next += real(conj(tmp[site][a]-source[site][a])*(tmp[site][a]-source[site][a]));
				}
			}
			reduceAllSum(norm_next);

			if (isOutputProcess()) std::cout << "LandauPropagator::Final norm for the ghost propagator after " << step << " steps: " << norm_next << std::endl;

			return true;
		}
		else {
			//std::cout << "ConjugateGradient::Residual norm at step " << step << ": " << norm_next << std::endl;
		}

		real_t beta = static_cast<real_t>(norm_next/norm);

#pragma omp parallel for
		for (int site = 0; site < source.completesize; ++site) {
			p[site] = r[site] + beta * p[site];
		}
		//p.updateHalo();
	}

	if (isOutputProcess()) std::cout << "LandauPropagator::Failure in finding convergence, last error: " << norm_next << std::endl;

	long_real_t result_re = 0., result_im = 0.;			
#pragma omp parallel for reduction(+:result_re,result_im)
	for (int site = 0; site < r.localsize; ++site) {
		real_t phase = 	   Layout::globalIndexX(site)*momentum[0]
					+ Layout::globalIndexY(site)*momentum[1] 
					+ Layout::globalIndexZ(site)*momentum[2] 
					+ Layout::globalIndexT(site)*momentum[3] ;
		
		result_re += real(std::complex<real_t>(cos(phase),-sin(phase))*ghost[site][c]);
		result_im += imag(std::complex<real_t>(cos(phase),-sin(phase))*ghost[site][c]);
	}
	reduceAllSum(result_re);
	reduceAllSum(result_im);
	result = std::complex<real_t>(result_re,result_im);

	return false;

}

void LandauPropagator::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;//TODO: only vector operations?
	typedef extended_gauge_lattice_t LT;

	std::vector<real_t> maximal_momentum = environment.configurations.get< std::vector<real_t> >("LandauPropagator::maximal_momentum");
	real_t cut =  environment.configurations.get< real_t >("LandauPropagator::anisotropy_cut");

	real_t epsilon = environment.configurations.get< real_t >("LandauPropagator::epsilon");
	unsigned int max_steps = environment.configurations.get< int >("LandauPropagator::max_steps");

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

	if (isOutputProcess()) std::cout << "LandauPropagator::Deviation from the Landau gauge: " << this->deviation(environment.gaugeLinkConfiguration) << std::endl;
	if (isOutputProcess()) std::cout << "LandauPropagator::Analizing " << momenta.size() << " momenta configurations ..." << std::endl;
	

	extended_gauge_lattice_t Afield;
	this->getLieAlgebraField(Afield, environment.gaugeLinkConfiguration);


	LieGenerator<GaugeGroup> lieGenerators;


	std::complex<real_t> Azero[4][lieGenerators.numberGenerators()];

	std::vector< std::vector<real_t> > data;

	for (unsigned int i = 0; i < momenta.size(); ++i) {
		std::vector<real_t> momentum = momenta[i];

		//First we compute the ghost propagator
		std::complex<real_t> ghost_propagator;
		this->ghostPropagator(ghost_propagator, environment.gaugeLinkConfiguration, momentum, 0, epsilon, max_steps);


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
					//std::cout << "Giusto per " << mu << " " << c << ": " << Azero[mu][c] << std::endl;
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
			//std::cout << "Boh: " << resultAk[0][0] * Azero[1][1] *resultAmk[0][2] * p[1] / 6. << std::endl;
			//std::cout << "Boh: " << resultAk[0][0] << " " << resultAmk[0][0] << " " << Azero[1][1] << std::endl;
		}

		if (isOutputProcess()) std::cout << "                     Gh2(p,-p): " << ghost_propagator << std::endl  << std::endl;
		result.push_back(real(ghost_propagator));

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
				int nn = 1;
				real_t averaged_result = data[i][1];
				for (unsigned int  j = i + 1; j < data.size(); ++j) {
					if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
						++nn;
						averaged_result += data[j][1];
					}
				}
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

	


		output->push("ghost_propagator");

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
				real_t averaged_result = data[i][3];
				for (unsigned int  j = i + 1; j < data.size(); ++j) {
					if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
						++nn;
						averaged_result += data[j][3];
					}
				}
				output->push("ghost_propagator");
				output->write("ghost_propagator", data[i][0]);
				output->write("ghost_propagator", averaged_result/nn);
				output->pop("ghost_propagator");
			}
		}

		output->pop("ghost_propagator");
	}

}

void LandauPropagator::registerParameters(po::options_description& desc) {
	desc.add_options()
		("LandauPropagator::epsilon", po::value<real_t>()->default_value(0.00000000001), "Convergence for the computation of the ghost propagator")
		("LandauPropagator::anisotropy_cut", po::value<real_t>()->default_value(0.25), "Cut for the off-diagonal momenta")
		("LandauPropagator::max_steps", po::value<int>()->default_value(5000), "Maximum inversion steps for the computation of the ghost propagator")
		("LandauPropagator::maximal_momentum", po::value<std::string>()->default_value("{2,2,2,2}"), "Momentum for the measure of the vertex function (syntax: {px,py,pz,pt})")
		;
}

}

