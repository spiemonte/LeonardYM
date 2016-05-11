/*
 * ConjugateGradient.cpp
 *
 *  Created on: Sep 27, 2012
 *      Author: spiem_01
 */

#include "ConjugateGradient.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

ConjugateGradient::ConjugateGradient() : epsilon(0.00000000001), maxSteps(3000) { }

ConjugateGradient::~ConjugateGradient() { }

bool ConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& original_source, reduced_dirac_vector_t& original_solution, reduced_dirac_vector_t const* initial_guess) {
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution;
	if (initial_guess == 0) {
		solution = source;
	} else {
		solution = *initial_guess;
	}
	
	dirac->multiply(tmp,solution);

#pragma omp parallel for
	for (int site = 0; site < source.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			r[site][mu] = source[site][mu] - tmp[site][mu];
			p[site][mu] = r[site][mu];
		}
	}

	long_real_t norm = AlgebraUtils::squaredNorm(r);
	
	long_real_t norm_next = norm;

	for (unsigned int step = 0; step < maxSteps; ++step) {
		dirac->multiply(tmp,p);
		norm = norm_next;
		std::complex<real_t> alpha = static_cast< std::complex<real_t> >(norm/AlgebraUtils::dot(p,tmp));


#pragma omp parallel for
		for (int site = 0; site < source.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + alpha * p[site][mu];
				r[site][mu] = r[site][mu] - alpha * tmp[site][mu];
			}
		}
		//solution.updateHalo();
		//r.updateHalo();//TODO maybe not needed

		norm_next = AlgebraUtils::squaredNorm(r);
		if (norm_next < epsilon) {
			lastSteps = step;
			original_solution = solution;
			return true;
		}
		else {
			//std::cout << "ConjugateGradient::Residual norm at step " << step << ": " << norm_next << std::endl;
		}

		real_t beta = static_cast<real_t>(norm_next/norm);

#pragma omp parallel for
		for (int site = 0; site < source.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = r[site][mu] + beta * p[site][mu];
			}
		}
		//p.updateHalo();
	}

	lastSteps = maxSteps;
	original_solution = solution;
	if (isOutputProcess()) std::cout << "ConjugateGradient::Failure in finding convergence, last error: " << norm_next << std::endl;
	return false;
}

void ConjugateGradient::setPrecision(double _epsilon) {
	epsilon = _epsilon;
}

double ConjugateGradient::getPrecision() const {
	return epsilon;
}

void ConjugateGradient::setMaximumSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

unsigned int ConjugateGradient::getMaximumSteps() const {
	return maxSteps;
}

double ConjugateGradient::getLastError() const {
	return lastError;
}

unsigned int ConjugateGradient::getLastSteps() const {
	return lastSteps;
}

} /* namespace Update */
