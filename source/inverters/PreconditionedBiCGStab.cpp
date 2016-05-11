/*
 * BiConjugateGradient.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#include "PreconditionedBiCGStab.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

PreconditionedBiCGStab::PreconditionedBiCGStab(const PreconditionedBiCGStab& toCopy) : Solver(toCopy), biConjugateGradient(0) { }

PreconditionedBiCGStab::PreconditionedBiCGStab() : Solver(), biConjugateGradient(0) { }

PreconditionedBiCGStab::~PreconditionedBiCGStab() {
	if (biConjugateGradient != 0) delete biConjugateGradient;
}

#ifdef ENABLE_MPI

bool PreconditionedBiCGStab::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, extended_dirac_vector_t& original_solution, DiracOperator* preconditioner, extended_dirac_vector_t const* original_initial_guess) {
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

bool PreconditionedBiCGStab::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* preconditioner, reduced_dirac_vector_t const* initial_guess) {
	solution = source;
	AlgebraUtils::gamma5(solution);
	reduced_dirac_vector_t tmp;
	dirac->multiply(tmp, solution);
	AlgebraUtils::gamma5(tmp);

	DiracOperator* squareDiracOperator = DiracOperator::getSquare(dirac);
	squareDiracOperator->setGamma5(true);
	

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(this->getMaximumSteps());
		biConjugateGradient->setPrecision(this->getPrecision());
	}

	bool res = biConjugateGradient->solve(squareDiracOperator, tmp, solution, initial_guess);
	delete squareDiracOperator;
	lastSteps = biConjugateGradient->getLastSteps();
	return res;
}

} /* namespace Update */
