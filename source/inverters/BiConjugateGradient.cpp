/*
 * BiConjugateGradient.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#include "BiConjugateGradient.h"
#include "algebra_utils/AlgebraUtils.h"
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

bool BiConjugateGradient::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution, extended_dirac_vector_t const* original_initial_guess) {
	//First set the initial solution
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution = source;
	bool result;
	if (original_initial_guess != 0) {
		reduced_dirac_vector_t initial_guess = *original_initial_guess;
		result = solve(dirac, source, solution, &initial_guess);
	}
	else result = solve(dirac, source, solution);
	original_solution = solution;
	return result;
}

#endif

bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess) {
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > epsilon) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use p as temporary vector
	dirac->multiply(p,solution);
	
	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - p[site][mu];
			residual_hat[site][mu] = source[site][mu] + p[site][mu];
		}
	}

	//Set p to residual
	//real_t beta = AlgebraUtils::squaredNorm(residual);

	//Set the initial parameter of the program
	unsigned int step = 0;
	std::complex<real_t> rho_prev, alpha, omega;
	
	while (step < maxSteps) {
		//rho = (residual_hat,residual)
		std::complex<real_t> rho = static_cast< std::complex<real_t> >(AlgebraUtils::dot(residual_hat,residual));
		if (rho.real()*rho.real() + rho.imag()*rho.imag() == 0.) return false;
		
		if (step == 0) {
			p = residual;
		}
		else {
			std::complex<real_t> beta = (rho/rho_prev)*(alpha/omega);
			//p = residual + beta*(p - omega*nu)
#pragma omp parallel for
			for (int site = 0; site < solution.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
				}
			}
		}
		//p_tilde = P.p
		preconditioner->multiply(p_tilde,p);
		
		//nu = A.p_tilde
		dirac->multiply(nu,p_tilde);
		
		//alpha = rho/(residual_hat,nu)
		alpha = rho/static_cast< std::complex<real_t> >(AlgebraUtils::dot(residual_hat,nu));
		
		//s = residual + alpha * nu
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha*nu[site][mu];
			}
		}
		
		
		//s_tilde = P.s
		preconditioner->multiply(s_tilde,s);
		
		//t = A.s
		dirac->multiply(t,s);
		
		//omega = (t,s)/(t,t)
		omega = static_cast< std::complex<real_t> >(AlgebraUtils::dot(t,s)/AlgebraUtils::squaredNorm(t));
		
		//x = x + alpha*p_tilde + omega*s_tilde
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + alpha*p_tilde[site][mu] + omega*s_tilde[site][mu];
			}
		}
		
		//residual = s - omega*t
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				residual[site][mu] = s[site][mu] - omega*t[site][mu];
			}
		}
		
		//Check convergence
		long_real_t norm = AlgebraUtils::squaredNorm(residual);
		if (norm < epsilon && step > 5) {
			lastSteps = step;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - final error norm: " << norm << std::endl;
#endif
			return true;
		}
		else {
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - error norm: " << norm << std::endl;
		}
		
		rho_prev = rho;
		
		lastError = norm;
		++step;
	}

	lastSteps = maxSteps;
	if (isOutputProcess()) std::cout << "BiConjugateGradient::Failure in finding convergence after " << maxSteps << " cicles, last error: " << lastError << std::endl;
	
	return false;
}

bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, int l, reduced_dirac_vector_t const* initial_guess) {
	typedef reduced_dirac_vector_t::Layout Layout;
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > epsilon) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use p as temporary vector
	dirac->multiply(p,solution);
	
	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < Layout::completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - p[site][mu];
			residual_hat[site][mu] = source[site][mu] + p[site][mu];
		}
	}

	//Set the initial parameter of the program
	unsigned int step = 0;
	std::complex<real_t> rho0 = 1., alpha = 0., omega = 1.;

	reduced_dirac_vector_t u_hat[l+1], r_hat[l+1];
	
	AlgebraUtils::setToZero(u_hat[0]);
	r_hat[0] = residual;
	
	while (step < maxSteps) {
		rho0 = -omega*rho0;

		for (int j = 0; j < l; ++j) {
			std::complex<real_t> rho1 = static_cast< std::complex<real_t> >(AlgebraUtils::dot(r_hat[j],residual_hat));
			std::complex<real_t> beta = (rho1/rho0)*alpha;
			rho0 = rho1;
			for (int i = 0; i <= j; ++i) {
#pragma omp parallel for
				for (int site = 0; site < Layout::completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						u_hat[i][site][mu] = r_hat[i][site][mu] - beta*u_hat[i][site][mu];
					}
				}
			}
			
			dirac->multiply(u_hat[j+1],u_hat[j]);

			std::complex<real_t> gamma = static_cast< std::complex<real_t> >(AlgebraUtils::dot(u_hat[j+1],residual_hat));
			alpha = rho0/gamma;

			for (int i = 0; i <= j; ++i) {
#pragma omp parallel for
				for (int site = 0; site < Layout::completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						r_hat[i][site][mu] = r_hat[i][site][mu] - alpha*u_hat[i+1][site][mu];
					}
				}
			}

			dirac->multiply(r_hat[j+1],r_hat[j]);

#pragma omp parallel for
			for (int site = 0; site < Layout::completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = solution[site][mu] + alpha*u_hat[0][site][mu];
				}
			}
			
		}

		std::complex<real_t> tau[l+1][l+1];
		std::complex<real_t> sigma[l+1];
		std::complex<real_t> gamma[l+1];
		std::complex<real_t> gamma_p[l+1];
		std::complex<real_t> gamma_pp[l+1];

		for (int j = 1; j < l + 1; ++j) {
			for (int i = 1; i < j; ++i) {
				tau[i][j] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(r_hat[j],r_hat[i]))/sigma[i];

#pragma omp parallel for
				for (int site = 0; site < Layout::completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						r_hat[j][site][mu] = r_hat[j][site][mu] - tau[i][j]*r_hat[i][site][mu];
					}
				}
			}
			
			sigma[j] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(r_hat[j],r_hat[j]));

			gamma_p[j] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(r_hat[0],r_hat[j]))/sigma[j];
		}

		gamma[l] = gamma_p[l];
		omega = gamma[l];

		for (int j = l-1; j > 0; --j) {
			std::complex<real_t> tmp = 0.;
			for (int i = j+1; i <= l; ++i) {
				tmp += tau[j][i]*gamma[i];
			}

			gamma[j] = gamma_p[j] - tmp;
		}

		for (int j = 1; j < l; ++j) {
			std::complex<real_t> tmp = 0.;
			for (int i = j+1; i < l; ++i) {
				tmp += tau[j][i]*gamma[i+1];
			}

			gamma_pp[j] = gamma[j+1] + tmp;
		}


#pragma omp parallel for
		for (int site = 0; site < Layout::completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + gamma[1]*r_hat[0][site][mu];
			}
		}

#pragma omp parallel for
		for (int site = 0; site < Layout::completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				r_hat[0][site][mu] = r_hat[0][site][mu] - gamma_p[l]*r_hat[l][site][mu];
			}
		}

#pragma omp parallel for
		for (int site = 0; site < Layout::completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				u_hat[0][site][mu] = u_hat[0][site][mu] - gamma[l]*u_hat[l][site][mu];
			}
		}

		for (int j = 1; j < l; ++j) {
#pragma omp parallel for
			for (int site = 0; site < Layout::completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					u_hat[0][site][mu] = u_hat[0][site][mu] - gamma[j]*u_hat[j][site][mu];
				}
			}
			
#pragma omp parallel for
			for (int site = 0; site < Layout::completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = solution[site][mu] + gamma_pp[j]*r_hat[j][site][mu];
				}
			}

#pragma omp parallel for
			for (int site = 0; site < Layout::completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					r_hat[0][site][mu] = r_hat[0][site][mu] - gamma_p[j]*r_hat[j][site][mu];
				}
			}
		}
		
		//Check convergence
		long_real_t norm = AlgebraUtils::squaredNorm(r_hat[0]);
		if (norm < epsilon) {
			lastSteps = step;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - final error norm: " << real(norm) << std::endl;
#endif
			return true;
		}
		else {
			std::cout << "BiCgStab(l)::Norm at step " << step << ": " << norm << std::endl;
		}
		
		lastError = norm;
		++step;
	}

	lastSteps = maxSteps;
	if (isOutputProcess()) std::cout << "BiConjugateGradient::Failure in finding convergence after " << maxSteps << " cicles, last error: " << lastError << std::endl;
	
	return false;
}

/*bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess) {
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > epsilon) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use p as temporary vector
	dirac->multiply(p,solution);
	
	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - p[site][mu];
			residual_hat[site][mu] = source[site][mu] + p[site][mu];
		}
	}

	//Set p to residual
	p = residual;

	//Set the initial parameter of the program
	unsigned int step = 0;

	while (step < maxSteps) {
		//p_tilde = P.p
		//preconditioner->multiply(p_tilde,p);
		p_tilde = p;
		
		//nu = A.p_tilde
		dirac->multiply(nu,p_tilde);
		//std::cout << AlgebraUtils::squaredNorm(nu) << std::endl;
		//alpha = (residual,residual_hat)/(nu,residual_hat)
		std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(residual,residual_hat)/AlgebraUtils::dot(nu,residual_hat));
		
		//s = residual - alpha*nu
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha*nu[site][mu];
			}
		}
		
		//s_tilde = P.s
		//preconditioner->multiply(s_tilde,s);
		s_tilde = s;
		
		//t = A.s_tilde
		dirac->multiply(t,s_tilde);
		
		//omega = (t,s)/(t,t)
		std::complex<real_t> omega = static_cast< std::complex<real_t> >(AlgebraUtils::dot(t,s)/AlgebraUtils::squaredNorm(t));
		
		//x = x + alpha*p_tilde + omega*s_tilde
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + alpha*p_tilde[site][mu] + omega*s_tilde[site][mu];
			}
		}
		
		//gamma = (residual,residual_hat)
		std::complex<real_t> gamma = static_cast< std::complex<real_t> >(AlgebraUtils::dot(residual,residual_hat));
		
		//r = s - omega * t
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				residual[site][mu] = s[site][mu] - omega*t[site][mu];
			}
		}
		
		long_real_t norm = AlgebraUtils::squaredNorm(residual);
		if (norm < epsilon) {
			lastSteps = step;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "BiCGStab steps: " << step << " - final error norm: " << real(norm) << std::endl;
#endif
			return true;
		}
		else {
			std::cout << "Norm at step " << step << " : " << norm << std::endl;
		}
		
		//beta = (residual,residual_hat)*alpha/(omega*gamma)
		std::complex<real_t> beta = static_cast< std::complex<real_t> >(AlgebraUtils::dot(residual,residual_hat))*alpha/(omega*gamma);
		
		//p = residual + beta*(p - omega*nu)
#pragma omp parallel for
		for (int site = 0; site < solution.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
			}
		}

		lastError = norm;
		++step;
	}

	lastSteps = maxSteps;
	if (isOutputProcess()) std::cout << "BiConjugateGradient::Failure in finding convergence after " << maxSteps << " cicles, last error: " << lastError << std::endl;
	
	return false;
}*/

bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess) {
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > epsilon) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use p as temporary vector
	dirac->multiply(p,solution);

	//Set the initial residual to source-A.solution and residual_hat accordingly
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - p[site][mu];
			residual_hat[site][mu] = source[site][mu] + p[site][mu];
		}
	}

	//residual.updateHalo();
	//residual_hat.updateHalo();

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
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
			}
		}
		//p.updateHalo();

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
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha *(nu[site][mu]);
			}
		}
		//s.updateHalo();

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
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = source[site][mu];
				}
			}
			//solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] += alpha*(p[site][mu]) + omega*(s[site][mu]);
			}
		}
		//solution.updateHalo();

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
		//residual.updateHalo();//TODO maybe not needed
#pragma omp parallel for
		for (int site = solution.localsize; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				residual[site][mu] = s[site][mu] - omega*(t[site][mu]);
			}
		}


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


bool BiConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, const std::complex<real_t>& shift, reduced_dirac_vector_t const* initial_guess) {
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > epsilon) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use p as temporary vector
	dirac->multiplyAdd(p,solution,solution,shift);

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
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*(p[site][mu] - omega*nu[site][mu]);
			}
		}
		//p.updateHalo();

		//nu = A.p[[k]]
		dirac->multiplyAdd(nu,p,p,shift);

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
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				s[site][mu] = residual[site][mu] - alpha *(nu[site][mu]);
			}
		}
		//s.updateHalo();

		//t = A.s;
		dirac->multiplyAdd(t,s,s,shift);

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
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = source[site][mu];
				}
			}
			//solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] += alpha*(p[site][mu]) + omega*(s[site][mu]);
			}
		}
		//solution.updateHalo();

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
