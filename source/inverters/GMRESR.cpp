/*
 * BiConjugateGradient.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#include "GMRESR.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

GMRESR::GMRESR() : Solver("GMRESR"), conjugateSpaceDimension(13) {
	c = new reduced_dirac_vector_t[conjugateSpaceDimension];
	u = new reduced_dirac_vector_t[conjugateSpaceDimension];
}

GMRESR::GMRESR(const GMRESR& toCopy) : Solver(toCopy), conjugateSpaceDimension(toCopy. conjugateSpaceDimension) {
	c = new reduced_dirac_vector_t[conjugateSpaceDimension];
	u = new reduced_dirac_vector_t[conjugateSpaceDimension];
}

GMRESR::~GMRESR() {
	delete[] c;
	delete[] u;
}

#ifdef ENABLE_MPI

bool GMRESR::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution, DiracOperator* preconditioner, extended_dirac_vector_t const* original_initial_guess) {
	//First set the initial solution
	reduced_dirac_vector_t source = original_source;
	reduced_dirac_vector_t solution = source;
	bool result;
	if (original_initial_guess != 0) {
		reduced_dirac_vector_t initial_guess = *original_initial_guess;
		result = solve(dirac, source, solution, preconditioner, &initial_guess);
	}
	else result = solve(dirac, source, solution, preconditioner, 0);
	original_solution = solution;
	return result;
}

#endif

bool GMRESR::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess) {
	//First set the initial solution
	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > precision) {
			solution = source;
		}
		else {
			AlgebraUtils::generateRandomVector(solution);
		}
	} else {
		solution = *initial_guess;
	}
	
	//Use r as temporary vector
	dirac->multiply(r,solution);
	
	//Set the initial residual to source-A.solution
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			r[site][mu] = source[site][mu] - r[site][mu];
		}
	}

	for (unsigned int k = 0; k < maxSteps; ++k) {
		if (preconditioner != 0) preconditioner->multiply(u[(k % conjugateSpaceDimension)],r);
		else u[(k % conjugateSpaceDimension)] = r;

		dirac->multiply(c[(k % conjugateSpaceDimension)],u[(k % conjugateSpaceDimension)]);
			
		for (unsigned int i = 0; i < (k % conjugateSpaceDimension); ++i) {
			std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(c[i],c[(k % conjugateSpaceDimension)]));
#pragma omp parallel for
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					c[(k % conjugateSpaceDimension)][site][mu] = c[(k % conjugateSpaceDimension)][site][mu] - alpha*c[i][site][mu];
					u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu] - alpha*u[i][site][mu];
				}
			}		
		}
			
		real_t norm = sqrt(AlgebraUtils::squaredNorm(c[(k % conjugateSpaceDimension)]));
#pragma omp parallel for
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				c[(k % conjugateSpaceDimension)][site][mu] = c[(k % conjugateSpaceDimension)][site][mu]/norm;
				u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu]/norm;
			}
		}
			
		std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(c[(k % conjugateSpaceDimension)],r));
#pragma omp parallel for
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + alpha*u[(k % conjugateSpaceDimension)][site][mu];
				r[site][mu] = r[site][mu] - alpha*c[(k % conjugateSpaceDimension)][site][mu];
			}
		}

		long_real_t error = AlgebraUtils::squaredNorm(r);
		if (error < precision) {
			lastSteps = k;
			return true;
		}
		else {
			//if (isOutputProcess()) std::cout << "GMRESR::Residual norm at step " << k << ": " << error << std::endl;
			lastError = error;
		}

	}

	lastSteps = maxSteps;
	if (isOutputProcess()) std::cout << "GMRESR::Failure in finding convergence after " << maxSteps << " cicles, last error: " << lastError << std::endl;
	
	return false;
}

void GMRESR::setConjugateSpaceDimension(unsigned int _conjugateSpaceDimension) {
	conjugateSpaceDimension = _conjugateSpaceDimension;
}

unsigned int GMRESR::getConjugateSpaceDimension() const {
	return conjugateSpaceDimension;
}

} /* namespace Update */
