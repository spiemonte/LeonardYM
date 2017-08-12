#include "LandauGhostPropagator.h"
#include "ToString.h"
#include "io/GlobalOutput.h"
#include "utils/LieGenerators.h"

namespace Update {

LandauGhostPropagator::LandauGhostPropagator() : LandauGaugeFixing() { }

LandauGhostPropagator::LandauGhostPropagator(const LandauGhostPropagator& toCopy) : LandauGaugeFixing(toCopy) { }

LandauGhostPropagator::~LandauGhostPropagator() { }

void LandauGhostPropagator::ghostMatrix(extended_adjoint_color_vector_t& output, const extended_adjoint_color_vector_t& input, const extended_adjoint_lattice_t& A, const extended_adjoint_lattice_t& B, const extended_adjoint_lattice_t& C) const {
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

bool LandauGhostPropagator::ghostPropagatorCG(std::complex<real_t>& result, const extended_gauge_lattice_t& lattice, const std::vector<real_t>& momentum, int c, real_t epsilon, unsigned int max_steps) const {
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
	source.updateHalo();

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

			if (isOutputProcess()) std::cout << "LandauGhostPropagator::Final norm for the ghost propagator after " << step << " steps: " << norm_next << std::endl;

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

	if (isOutputProcess()) std::cout << "LandauGhostPropagator::Failure in finding convergence, last error: " << norm_next << std::endl;

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

bool LandauGhostPropagator::ghostPropagatorBiCGStab(std::complex<real_t>& result, const extended_gauge_lattice_t& lattice, const std::vector<real_t>& momentum, int c, real_t epsilon, unsigned int max_steps) const {
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
	extended_adjoint_color_vector_t source, tmp, residual, p, s, t, nu, residual_hat;
#pragma omp parallel for
	for (int site = 0; site < source.localsize; ++site) {
		real_t phase = 	   Layout::globalIndexX(site)*momentum[0]
					+ Layout::globalIndexY(site)*momentum[1] 
					+ Layout::globalIndexZ(site)*momentum[2] 
					+ Layout::globalIndexT(site)*momentum[3] ;
		set_to_zero(source[site]);
		source[site][c] = std::complex<real_t>(cos(phase),sin(phase));
	}
	source.updateHalo();

	ghost = source;
	
	//Use p as temporary vector
	this->ghostMatrix(p,ghost,A,B,C);

	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < ghost.completesize; ++site) {
		residual[site] = source[site] - p[site];
		residual_hat[site] = source[site] + p[site];
	}

	//Set nu and p to zero
#pragma omp parallel for
	for (int site = 0; site < ghost.completesize; ++site) {
		set_to_zero(p[site]);
		set_to_zero(nu[site]);
	}

	//Set the initial parameter of the program
	std::complex<real_t> alpha = 1., omega = 1.;
	std::complex<long_real_t> rho = 1.;
	unsigned int step = 0;
	real_t minimal_norm = 10000.;

	while (step < max_steps) {
		//rho[k] = rhat.r[k-1]
		long_real_t rho_next_re = 0.;
		long_real_t rho_next_im = 0.;
#pragma omp parallel for reduction(+:rho_next_re, rho_next_im)
		for (int site = 0; site < ghost.localsize; ++site) {
			std::complex<real_t> partial = vector_dot(residual_hat[site],residual[site]);
			rho_next_re += real(partial);
			rho_next_im += imag(partial);
		}
		reduceAllSum(rho_next_re);
		reduceAllSum(rho_next_im);

		std::complex<long_real_t> rho_next(rho_next_re,rho_next_im);

		if (norm(rho_next) == 0.) {
			if (isOutputProcess()) std::cout << "LandauGhostPropagator::Fatal error in norm " << rho_next << " at step " << step << " in BiCGStab!"<< std::endl;
			return false;//TODO
		}

		std::complex<real_t> beta = static_cast< std::complex<real_t> >((rho_next/rho))*(alpha/omega);
		//p = r[[k - 1]] + beta*(p[[k - 1]] - omega[[k - 1]]*nu[[k - 1]])
#pragma omp parallel for
		for (int site = 0; site < ghost.completesize; ++site) {
			p[site] = residual[site] + beta*(p[site] - omega*nu[site]);
		}
		//p.updateHalo();

		//nu = A.p[[k]]
		this->ghostMatrix(nu,p,A,B,C);

		//alpha = rho[[k]]/(rhat[[1]].nu[[k]]);
		long_real_t alphatmp_re = 0.;
		long_real_t alphatmp_im = 0.;
#pragma omp parallel for reduction(+:alphatmp_re, alphatmp_im)
		for (int site = 0; site < ghost.localsize; ++site) {
			std::complex<real_t> partial = vector_dot(residual_hat[site],nu[site]);
			alphatmp_re += real(partial);
			alphatmp_im += imag(partial);
		}
		reduceAllSum(alphatmp_re);
		reduceAllSum(alphatmp_im);

		std::complex<long_real_t> alphatmp(alphatmp_re,alphatmp_im);
		alpha = static_cast< std::complex<real_t> >(rho_next/alphatmp);

		//s = r[[k - 1]] - alpha*nu[[k]]
#pragma omp parallel for
		for (int site = 0; site < ghost.completesize; ++site) {
			s[site] = residual[site] - alpha *(nu[site]);
		}
		//s.updateHalo();

		//t = A.s;
		this->ghostMatrix(t,s,A,B,C);

		//omega = (t.s)/(t.t)
		long_real_t tmp1_re = 0., tmp1_im = 0., tmp2_re = 0., tmp2_im = 0.;
#pragma omp parallel for reduction(+:tmp1_re, tmp1_im, tmp2_re, tmp2_im)
		for (int site = 0; site < ghost.localsize; ++site) {
			std::complex<real_t> partial1 = vector_dot(t[site],s[site]);
			tmp1_re += real(partial1);
			tmp1_im += imag(partial1);
			std::complex<real_t> partial2 = vector_dot(t[site],t[site]);
			tmp2_re += real(partial2);
			tmp2_im += imag(partial2);
		}
		reduceAllSum(tmp1_re);
		reduceAllSum(tmp1_im);
		reduceAllSum(tmp2_re);
		reduceAllSum(tmp2_im);

		std::complex<long_real_t> tmp1(tmp1_re, tmp1_im), tmp2(tmp2_re, tmp2_im);
		omega = static_cast< std::complex<real_t> >(tmp1/tmp2);

		if (real(tmp2) == 0) {
#pragma omp parallel for
			for (int site = 0; site < ghost.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					ghost[site][mu] = source[site][mu];
				}
			}
			//solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
		for (int site = 0; site < ghost.completesize; ++site) {
			ghost[site] += alpha*(p[site]) + omega*(s[site]);
		}
		//solution.updateHalo();

		//residual[[k]] = s - omega[[k]]*t
		//norm = residual[[k]].residual[[k]]
		long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
		for (int site = 0; site < ghost.localsize; ++site) {
			residual[site] = s[site] - omega*(t[site]);
			norm += real(vector_dot(residual[site],residual[site]));
		}
		reduceAllSum(norm);
		residual.updateHalo();//TODO maybe not needed

		if (norm < epsilon) {
			long_real_t result_re = 0., result_im = 0.;			
#pragma omp parallel for reduction(+:result_re,result_im)
			for (int site = 0; site < ghost.localsize; ++site) {
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
			long_real_t norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
			for (int site = 0; site < ghost.localsize; ++site) {
				for (unsigned int a = 0; a < lieGenerators.numberGenerators(); ++a) {
					norm_next += real(conj(tmp[site][a]-source[site][a])*(tmp[site][a]-source[site][a]));
				}
			}
			reduceAllSum(norm_next);

			if (isOutputProcess()) std::cout << "LandauGhostPropagator::Final norm for the ghost propagator after " << step << " steps: " << norm_next << std::endl;

			return true;
		}
		else if (norm < minimal_norm) {
			minimal_norm = norm;
			
			long_real_t result_re = 0., result_im = 0.;			
#pragma omp parallel for reduction(+:result_re,result_im)
			for (int site = 0; site < ghost.localsize; ++site) {
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
		}
//#ifdef BICGLOG
		//else if (isOutputProcess()) std::cout << "Error at step " << step << ": " << norm << std::endl;
//#endif

		rho = rho_next;

		++step;
	}

	if (isOutputProcess()) std::cout << "LandauGhostPropagator::Failure in finding convergence after " << max_steps << " bicg cicles, minimal norm found is " << minimal_norm << std::endl;
	
	
	
	return false;

}

void LandauGhostPropagator::execute(environment_t& environment) {
	typedef extended_gauge_lattice_t::Layout Layout;//TODO: only vector operations?
	typedef extended_gauge_lattice_t LT;

	std::vector<real_t> maximal_momentum = environment.configurations.get< std::vector<real_t> >("LandauGhostPropagator::maximal_momentum");
	real_t cut =  environment.configurations.get< real_t >("LandauGhostPropagator::anisotropy_cut");

	real_t epsilon = environment.configurations.get< real_t >("LandauGhostPropagator::epsilon");
	unsigned int max_steps = environment.configurations.get< int >("LandauGhostPropagator::max_steps");

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

	if (isOutputProcess()) std::cout << "LandauGhostPropagator::Deviation from the Landau gauge: " << convergence << std::endl;
	if (isOutputProcess()) std::cout << "LandauGhostPropagator::Analizing " << momenta.size() << " momenta configurations ..." << std::endl;
	

	extended_gauge_lattice_t Afield;
	this->getLieAlgebraField(Afield, environment.gaugeLinkConfiguration);


	LieGenerator<GaugeGroup> lieGenerators;


	std::complex<real_t> Azero[4][lieGenerators.numberGenerators()];

	std::vector< std::vector<real_t> > data;

	for (unsigned int i = 0; i < momenta.size(); ++i) {
		std::vector<real_t> momentum = momenta[i];

		//Here we compute the ghost propagator
		std::complex<real_t> ghost_propagator = 0.;
		if (i != 0) this->ghostPropagatorBiCGStab(ghost_propagator, environment.gaugeLinkConfiguration, momentum, 0, epsilon, max_steps);

		//Physical momentum
		std::vector<real_t> p(4);
		p[0] = 2.*sin(momentum[0]/2.);
		p[1] = 2.*sin(momentum[1]/2.);
		p[2] = 2.*sin(momentum[2]/2.);
		p[3] = 2.*sin(momentum[3]/2.);

		real_t modulus = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);

		std::vector<real_t> result;
		result.push_back(modulus);

		if (isOutputProcess()) std::cout << "                     Gh2(p,-p): " << ghost_propagator << std::endl  << std::endl;
		result.push_back(real(ghost_propagator));

		data.push_back(result);
	}

	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

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
				real_t averaged_result = data[i][1];
				for (unsigned int  j = i + 1; j < data.size(); ++j) {
					if (fabs(data[i][0] - data[j][0]) < 0.0000000000000001) {
						++nn;
						averaged_result += data[j][1];
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

void LandauGhostPropagator::registerParameters(po::options_description& desc) {
	desc.add_options()
		("LandauGhostPropagator::epsilon", po::value<real_t>()->default_value(0.00000000001), "Convergence for the computation of the ghost propagator")
		("LandauGhostPropagator::anisotropy_cut", po::value<real_t>()->default_value(0.25), "Cut for the off-diagonal momenta")
		("LandauGhostPropagator::max_steps", po::value<int>()->default_value(5000), "Maximum inversion steps for the computation of the ghost propagator")
		("LandauGhostPropagator::maximal_momentum", po::value<std::string>()->default_value("{2,2,2,2}"), "Momentum for the measure of the vertex function (syntax: {px,py,pz,pt})")
		;
}

}

