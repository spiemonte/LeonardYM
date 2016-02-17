#include "MultiGridConjugateGradient.h"

namespace Update {

MultiGridConjugateGradientSolver::MultiGridConjugateGradientSolver() : precision(0.000001), maxSteps(300) { }

bool MultiGridConjugateGradientSolver::solve(MultiGridOperator* multiGridOperator, const multigrid_vector_t& source_hat, multigrid_vector_t& solution_hat) {
	solution_hat = source_hat;
	multiGridOperator->multiply(tmp_hat,solution_hat);

#pragma omp parallel for
	for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
		r_hat[m] = source_hat[m] - tmp_hat[m];
		p_hat[m] = r_hat[m];
	}

	long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
	for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
		norm += real(conj(r_hat[m])*r_hat[m]);
	}
	reduceAllSum(norm);
	long_real_t norm_next = norm;

	for (int innerStep = 0; innerStep < maxSteps; ++innerStep) {
		multiGridOperator->multiply(tmp_hat,p_hat);
		norm = norm_next;
		long_real_t gammaRe = 0;
		long_real_t gammaIm = 0;
#pragma omp parallel for reduction(+:gammaRe,gammaIm)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			std::complex<real_t> result = conj(p_hat[m])*tmp_hat[m];
			gammaRe += real(result);
			gammaIm += imag(result);
		}
		reduceAllSum(gammaRe);
		reduceAllSum(gammaIm);
		std::complex<real_t> alpha = static_cast<real_t>(norm)/std::complex<real_t>(gammaRe,gammaIm);


#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			solution_hat[m] = solution_hat[m] + alpha * p_hat[m];
			r_hat[m] = r_hat[m] - alpha * tmp_hat[m];
		}

		norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			norm_next += real(conj(r_hat[m])*r_hat[m]);
		}
		//Check for convergence
		if (norm_next < precision && innerStep > 5) {
			if (isOutputProcess()) std::cout << "Inner convergence in " << innerStep << " steps."<< std::endl;
			break;
		} else if (innerStep == maxSteps - 1) {
			if (isOutputProcess()) std::cout << "Failure in finding convergence in inner solver, last error " << norm_next << std::endl;
			return false;
		}
		//std::cout << "Inner norm at step " << innerStep << ": " << norm_next << std::endl; 

		real_t beta = static_cast<real_t>(norm_next/norm);

#pragma omp parallel for
		for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
			p_hat[m] = r_hat[m] + beta * p_hat[m];
		}
	}

	return true;

}

void MultiGridConjugateGradientSolver::setPrecision(const real_t& _precision) {
	precision = _precision;
}

void MultiGridConjugateGradientSolver::setMaximumSteps(int _maxSteps) {
	maxSteps = _maxSteps;
}

}
