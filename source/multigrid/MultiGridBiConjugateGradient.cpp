#include "MultiGridBiConjugateGradient.h"

namespace Update {

MultiGridBiConjugateGradientSolver::MultiGridBiConjugateGradientSolver() : precision(0.000001), maxSteps(300) { }

bool MultiGridBiConjugateGradientSolver::solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source, multigrid_vector_t& solution) {
	solution = source;
	//Use p as temporary vector
	multiGridOperator->multiply(p,solution);
		
#pragma omp parallel for
	for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
		//Set the initial residual to source-A.solution and residual_hat accordingly
		residual[m] = source[m] - p[m];
		residual_hat[m] = source[m] + p[m];
		//Set then nu and p to zero
		p[m] = 0;
		nu[m] = 0;
	}

	//Set the initial parameter of the program
	std::complex<real_t> alpha = 1., omega = 1.;
	std::complex<long_real_t> rho = 1.;
	unsigned int step = 0;
	long_real_t norm_r = 0.;

	while (step < maxSteps) {
		//rho[k] = rhat.r[k-1]
		long_real_t rho_next_re = 0.;
		long_real_t rho_next_im = 0.;
#pragma omp parallel for reduction(+:rho_next_re, rho_next_im)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			complex partial = conj(residual_hat[m])*residual[m];
			rho_next_re += real(partial);
			rho_next_im += imag(partial);
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
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			p[m] = residual[m] + beta*(p[m] - omega*nu[m]);
		}

		//nu = A.p[[k]]
		multiGridOperator->multiply(nu,p);

		//alpha = rho[[k]]/(rhat[[1]].nu[[k]]);
		long_real_t alphatmp_re = 0.;
		long_real_t alphatmp_im = 0.;
#pragma omp parallel for reduction(+:alphatmp_re, alphatmp_im)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			complex partial = conj(residual_hat[m])*nu[m];
			alphatmp_re += real(partial);
			alphatmp_im += imag(partial);
		}
		reduceAllSum(alphatmp_re);
		reduceAllSum(alphatmp_im);

		std::complex<long_real_t> alphatmp(alphatmp_re,alphatmp_im);
		alpha = static_cast< std::complex<real_t> >(rho_next/alphatmp);

		//s = r[[k - 1]] - alpha*nu[[k]]
#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			s[m] = residual[m] - alpha *(nu[m]);
		}
		//s.updateHalo();

		//t = A.s;
		multiGridOperator->multiply(t,s);

		//omega = (t.s)/(t.t)
		long_real_t tmp1_re = 0., tmp1_im = 0., tmp2_re = 0., tmp2_im = 0.;
#pragma omp parallel for reduction(+:tmp1_re, tmp1_im, tmp2_re, tmp2_im)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			complex partial1 = conj(t[m])*s[m];
			tmp1_re += real(partial1);
			tmp1_im += imag(partial1);
			complex partial2 = conj(t[m])*t[m];
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
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				solution[m] = source[m];
			}
			//solution.updateHalo();
			return true;//TODO, identity only?
		}

		//solution[[k]] = solution[[k - 1]] + alpha*p[[k]] + omega[[k]]*s
#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			solution[m] += alpha*(p[m]) + omega*(s[m]);
		}
		//solution.updateHalo();

		//residual[[k]] = s - omega[[k]]*t
		//norm = residual[[k]].residual[[k]]
		norm_r = 0.;
#pragma omp parallel for reduction(+:norm_r)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			residual[m] = s[m] - omega*(t[m]);
			norm_r += real(conj(residual[m])*residual[m]);
		}
		reduceAllSum(norm_r);


		if (norm_r < precision && step > 4) {
			//lastSteps = step;
#ifdef BICGLOG
			if (isOutputProcess()) std::cout << "MultiGridBiConjugateGradientSolver::Convergence in " << step << " - final error norm: " << norm_r << std::endl;
#endif
			return true;
		}
		//#ifdef BICGLOG
		//		else if (isOutputProcess()) std::cout << "Error at step " << step << ": " << real(norm) << std::endl;
		//#endif

		rho = rho_next;

		++step;
	}

#ifdef BICGLOG
	if (isOutputProcess()) std::cout << "MultiGridBiConjugateGradientSolver::Failure in finding convergence in " << maxSteps << " - final error norm: " << norm_r << std::endl;
#endif

	return false;
}

void MultiGridBiConjugateGradientSolver::setPrecision(const real_t& _precision) {
	precision = _precision;
}

void MultiGridBiConjugateGradientSolver::setMaximumSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

}


