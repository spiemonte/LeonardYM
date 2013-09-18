/*
 * BiConjugateGradient.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#include "BiConjugateGradient.h"
#include "AlgebraUtils.h"
//#define BICGLOG

namespace Update {

BiConjugateGradient::BiConjugateGradient() : epsilon(0.00000000001), maxSteps(3000) { }

BiConjugateGradient::~BiConjugateGradient() { }

/*bool BiConjugateGradient::solve(DiracOperator* dirac, const dirac_vector_t& source, dirac_vector_t& solution) {
	std::cout << "Norma soluzione: " << AlgebraUtils::squaredNorm(source) << std::endl;
	//First set the initial solution, nu and p to zero
#ifdef MULTITHREADING
	#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < solution.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			solution[site][mu].zeros();
			p[site][mu].zeros();
			nu[site][mu].zeros();
		}
	}

	solution.updateHalo();
	p.updateHalo();
	nu.updateHalo();

	//Set the initial residual to source-A.solution (that is to source) and residual_hat accordingly
#ifdef MULTITHREADING
	#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < solution.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu];
			residual_hat[site][mu] = residual[site][mu];
		}
	}

	residual.updateHalo();
	residual_hat.updateHalo();

	//Set the initial parameter of the program
	std::complex<real_t> alpha = 1., omega = 1.;
	std::complex<long_real_t> rho = 1.;
	unsigned int step = 0;

	while (step < maxSteps) {
		//rho[k] = rhat.r[k-1]
		long_real_t rho_next_re = 0.;
		long_real_t rho_next_im = 0.;
#ifdef MULTITHREADING
		#pragma omp parallel for reduction(+:rho_next_re, rho_next_im)
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(residual_hat[site][mu],residual[site][mu]);
				rho_next_re += real(partial);
				rho_next_im += imag(partial);
			}
		}
#ifdef ENABLE_MPI
		::MpiExt::reduceAllSum(rho_next_re);
		::MpiExt::reduceAllSum(rho_next_im);
#endif
		std::complex<long_real_t> rho_next(rho_next_re,rho_next_im);

		if (norm(rho_next) == 0.) {
			std::cout << resisula
			if (isOutputProcess()) std::cout << "BiConjugateGradient::Fatal error in norm " << rho_next << " at step " << step << " " << AlgebraUtils::dot(residual_hat,residual) << " " << AlgebraUtils::squaredNorm(residual_hat) << " " << AlgebraUtils::squaredNorm(residual) << std::endl;
			return false;//TODO
		}

		std::complex<real_t> beta = static_cast< std::complex<real_t> >((rho_next/rho))*(alpha/omega);
		//p = r[[k - 1]] + beta*(p[[k - 1]] - omega[[k - 1]]*nu[[k - 1]])
#ifdef MULTITHREADING
		#pragma omp parallel for
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
			}
		}
		p.updateHalo();

		//nu = A.p[[k]]
		dirac->multiply(nu,p);

		//alpha = rho[[k]]/(rhat[[1]].nu[[k]]);
		long_real_t alphatmp_re = 0.;
		long_real_t alphatmp_im = 0.;
#ifdef MULTITHREADING
		#pragma omp parallel for reduction(+:alphatmp_re, alphatmp_im)
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(residual_hat[site][mu],nu[site][mu]);
				alphatmp_re += real(partial);
				alphatmp_im += imag(partial);
			}
		}
#ifdef ENABLE_MPI
		::MpiExt::reduceAllSum(alphatmp_re);
		::MpiExt::reduceAllSum(alphatmp_im);
#endif
		std::complex<long_real_t> alphatmp(alphatmp_re,alphatmp_im);
		alpha = static_cast< std::complex<real_t> >(rho_next/alphatmp);

		//s = r[[k - 1]] - alpha*nu[[k]]
#ifdef MULTITHREADING
		#pragma omp parallel for
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha *(nu[site][mu]);
			}
		}
		s.updateHalo();

		//t = A.s;
		dirac->multiply(t,s);

		//omega = (t.s)/(t.t)
		long_real_t tmp1_re = 0., tmp1_im = 0., tmp2_re = 0., tmp2_im = 0.;
#ifdef MULTITHREADING
		#pragma omp parallel for reduction(+:tmp1_re, tmp1_im, tmp2_re, tmp2_im)
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial1 = vector_dot(t[site][mu],s[site][mu]);
				tmp1_re += real(partial1);
				tmp1_im += imag(partial1);
				complex partial2 = vector_dot(t[site][mu],t[site][mu]);
				tmp2_re += real(partial2);
				tmp2_im += imag(partial2);
			}
		}
#ifdef ENABLE_MPI
		::MpiExt::reduceAllSum(tmp1_re);
		::MpiExt::reduceAllSum(tmp1_im);
		::MpiExt::reduceAllSum(tmp2_re);
		::MpiExt::reduceAllSum(tmp2_im);
#endif
		std::complex<long_real_t> tmp1(tmp1_re, tmp1_im), tmp2(tmp2_re, tmp2_im);
		omega = static_cast< std::complex<real_t> >(tmp1/tmp2);

		if (real(tmp2) == 0) {
#ifdef MULTITHREADING
			#pragma omp parallel for
#endif
			for (unsigned int site = 0; site < solution.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = source[site][mu];
				}
			}
			solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#ifdef MULTITHREADING
		#pragma omp parallel for
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] += alpha*(p[site][mu]) + omega*(s[site][mu]);
			}
		}
		solution.updateHalo();

		//residual[[k]] = s - omega[[k]]*t
		//norm = residual[[k]].residual[[k]]
		long_real_t norm = 0.;
#ifdef MULTITHREADING
		#pragma omp parallel for reduction(+:norm)
#endif
		for (unsigned int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				residual[site][mu] = s[site][mu] - omega*(t[site][mu]);
				norm += real(vector_dot(residual[site][mu],residual[site][mu]));
			}
		}
#ifdef ENABLE_MPI
		::MpiExt::reduceAllSum(norm);
#endif
		residual.updateHalo();//TODO maybe not needed


		if (norm < epsilon) {
			lastSteps = step;
			lastError = norm;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - final error norm: " << real(norm) << std::endl;
#endif
			return true;
		}
//#ifdef BICGLOG
//		else if (isOutputProcess()) std::cout << "Error at step " << step << ": " << real(norm) << std::endl;
//#endif

		rho = rho_next;

		++step;
	}

	std::cout << "Failure in finding convergence after " << maxSteps << " cicles!" << std::endl;
	return false;
}*/

#ifdef ENABLE_MPI

bool BiConjugateGradient::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution) {
	//First set the initial solution
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution = source;
	bool result = solve(dirac, source, solution);
	original_solution = solution;
	return result;
}

#endif



bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution) {
	//First set the initial solution
	solution = source;
	
	//Use p as temporary vector
	dirac->multiply(p,solution);

	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < solution.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - p[site][mu];
			residual_hat[site][mu] = source[site][mu] + p[site][mu];
		}
	}

	residual.updateHalo();
	residual_hat.updateHalo();

	//Set nu and p to zero
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(p[site][mu]);
			set_to_zero(nu[site][mu]);
		}
	}

	//Set the initial parameter of the program
	std::complex<real_t> alpha = 1., omega = 1.;
	std::complex<long_real_t> rho = 1.;
	unsigned int step = 0;

	while (step < maxSteps) {
		//rho[k] = rhat.r[k-1]
		long_real_t rho_next_re = 0.;
		long_real_t rho_next_im = 0.;
#pragma omp parallel for reduction(+:rho_next_re, rho_next_im)
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(residual_hat[site][mu],residual[site][mu]);
				rho_next_re += real(partial);
				rho_next_im += imag(partial);
			}
		}
		reduceAllSum(rho_next_re);
		reduceAllSum(rho_next_im);

		std::complex<long_real_t> rho_next(rho_next_re,rho_next_im);

		if (norm(rho_next) == 0.) {
			if (isOutputProcess()) std::cout << "BiConjugateGradient::Fatal error in norm " << rho_next << " at step " << step << std::endl;
			return false;//TODO
		}

		std::complex<real_t> beta = static_cast< std::complex<real_t> >((rho_next/rho))*(alpha/omega);
		//p = r[[k - 1]] + beta*(p[[k - 1]] - omega[[k - 1]]*nu[[k - 1]])
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
			}
		}
		p.updateHalo();

		//nu = A.p[[k]]
		dirac->multiply(nu,p);

		//alpha = rho[[k]]/(rhat[[1]].nu[[k]]);
		long_real_t alphatmp_re = 0.;
		long_real_t alphatmp_im = 0.;
#pragma omp parallel for reduction(+:alphatmp_re, alphatmp_im)
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial = vector_dot(residual_hat[site][mu],nu[site][mu]);
				alphatmp_re += real(partial);
				alphatmp_im += imag(partial);
			}
		}
		reduceAllSum(alphatmp_re);
		reduceAllSum(alphatmp_im);

		std::complex<long_real_t> alphatmp(alphatmp_re,alphatmp_im);
		alpha = static_cast< std::complex<real_t> >(rho_next/alphatmp);

		//s = r[[k - 1]] - alpha*nu[[k]]
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha *(nu[site][mu]);
			}
		}
		s.updateHalo();

		//t = A.s;
		dirac->multiply(t,s);

		//omega = (t.s)/(t.t)
		long_real_t tmp1_re = 0., tmp1_im = 0., tmp2_re = 0., tmp2_im = 0.;
#pragma omp parallel for reduction(+:tmp1_re, tmp1_im, tmp2_re, tmp2_im)
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				complex partial1 = vector_dot(t[site][mu],s[site][mu]);
				tmp1_re += real(partial1);
				tmp1_im += imag(partial1);
				complex partial2 = vector_dot(t[site][mu],t[site][mu]);
				tmp2_re += real(partial2);
				tmp2_im += imag(partial2);
			}
		}
		reduceAllSum(tmp1_re);
		reduceAllSum(tmp1_im);
		reduceAllSum(tmp2_re);
		reduceAllSum(tmp2_im);

		std::complex<long_real_t> tmp1(tmp1_re, tmp1_im), tmp2(tmp2_re, tmp2_im);
		omega = static_cast< std::complex<real_t> >(tmp1/tmp2);

		if (real(tmp2) == 0) {
#pragma omp parallel for
			for (int site = 0; site < solution.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = source[site][mu];
				}
			}
			solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] += alpha*(p[site][mu]) + omega*(s[site][mu]);
			}
		}
		solution.updateHalo();

		//residual[[k]] = s - omega[[k]]*t
		//norm = residual[[k]].residual[[k]]
		long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				residual[site][mu] = s[site][mu] - omega*(t[site][mu]);
				norm += real(vector_dot(residual[site][mu],residual[site][mu]));
			}
		}
		reduceAllSum(norm);
		residual.updateHalo();//TODO maybe not needed


		if (norm < epsilon) {
			lastSteps = step;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - final error norm: " << real(norm) << std::endl;
#endif
			return true;
		}
//#ifdef BICGLOG
//		else if (isOutputProcess()) std::cout << "Error at step " << step << ": " << real(norm) << std::endl;
//#endif

		rho = rho_next;

		lastError = norm;
		++step;
	}

	lastSteps = maxSteps;
	if (isOutputProcess()) std::cout << "Failure in finding convergence after " << maxSteps << " cicles, last error: " << lastError << std::endl;
	
	return false;
}

void BiConjugateGradient::setPrecision(double _epsilon) {
	epsilon = _epsilon;
}

double BiConjugateGradient::getPrecision() const {
	return epsilon;
}

void BiConjugateGradient::setMaximumSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

unsigned int BiConjugateGradient::getMaximumSteps() const {
	return maxSteps;
}

double BiConjugateGradient::getLastError() const {
	return lastError;
}

unsigned int BiConjugateGradient::getLastSteps() const {
	return lastSteps;
}

} /* namespace Update */
