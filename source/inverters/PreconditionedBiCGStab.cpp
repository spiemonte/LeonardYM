/*
 * BiConjugateGradient.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: spiem_01
 */

#include "PreconditionedBiCGStab.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_operators/SquareEvenOddImprovedDiracWilsonOperator.h"

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
	typedef reduced_dirac_vector_t::Layout Layout;
	bool res = false;

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(this->getMaximumSteps());
		biConjugateGradient->setPrecision(this->getPrecision());
	}

	if (dynamic_cast<ImprovedDiracWilsonOperator*>(dirac) && (Layout::glob_x % 2 == 0) && (Layout::glob_y % 2 == 0) && (Layout::glob_z % 2 == 0) && (Layout::glob_t % 2 == 0)) {
		ImprovedDiracWilsonOperator* dirac_improved = dynamic_cast<ImprovedDiracWilsonOperator*>(dirac);
		
		SquareEvenOddImprovedDiracWilsonOperator* eo_sq = new SquareEvenOddImprovedDiracWilsonOperator();
		eo_sq->setKappa(dirac_improved->getKappa());
		eo_sq->setCSW(dirac_improved->getCSW());
		eo_sq->setLattice(*dirac->getLattice());

		EvenOddImprovedDiracWilsonOperator* eo = new EvenOddImprovedDiracWilsonOperator();
		eo->setKappa(dirac_improved->getKappa());
		eo->setCSW(dirac_improved->getCSW());
		eo->setLattice(*dirac->getLattice());

		/*reduced_dirac_vector_t prd1 = source, prd2;
		eo->multiply(prd2, source);
		std::cout << "Vediamo: " << AlgebraUtils::differenceNorm(source,  prd2) << " " <<  AlgebraUtils::dot(source, source)  << std::endl;
		prd2 = source;
		eo->multiplyEvenEvenInverse(prd2);
		eo->multiplyOddOdd(prd2, EVEN);
		std::cout << "Vediamo: " << AlgebraUtils::differenceNorm(source,  prd2) << " " <<  AlgebraUtils::dot(prd2, source)  << std::endl;
		std::cout << "1) " << prd2[1][0] << std::endl;
		std::cout << "2) " << prd2[1][1] << std::endl;
		std::cout << "3) " << prd2[1][2] << std::endl;
		std::cout << "4) " << prd2[1][3] << std::endl;

		std::cout << "1) " << prd2[30][0] << std::endl;
		std::cout << "2) " << prd2[30][1] << std::endl;
		std::cout << "3) " << prd2[30][2] << std::endl;
		std::cout << "4) " << prd2[30][3] << std::endl;
		exit(5);*/
		/*eo->multiply(prd1, tmp);
		std::cout << "Vediamo: " << AlgebraUtils::dot(prd1, source) << std::endl;*/
	
		reduced_dirac_vector_t even_inverse, odd_source = source, odd_source_eo, odd_inverse;
		//We construct the source for the odd
		eo->multiplyEvenEvenInverse(odd_source);
		eo->multiplyEvenOdd(odd_source, odd_source, ODD);
#pragma omp parallel for
		 for (int site = 0; site < Layout::completesize; ++site) {
			if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == 1) {
				for (unsigned int mu = 0; mu < 4; ++mu) odd_source[site][mu] = source[site][mu] - odd_source[site][mu];
			}
			else {
				for (unsigned int mu = 0; mu < 4; ++mu) odd_source[site][mu] = source[site][mu];
			}
			for (unsigned int mu = 2; mu < 4; ++mu) odd_source[site][mu] = -odd_source[site][mu];
		}
	
		//We take the inverse for the odd_source
		eo->multiply(odd_source_eo, odd_source);
		res = biConjugateGradient->solve(eo_sq, odd_source_eo, odd_inverse, initial_guess);	
	
		//Inverse for the even
		eo->multiplyEvenOdd(even_inverse, odd_inverse, EVEN);
#pragma omp parallel for
		for (int site = 0; site < Layout::completesize; ++site) {
			if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == 0) {
				for (unsigned int mu = 0; mu < 4; ++mu) even_inverse[site][mu] = source[site][mu] - even_inverse[site][mu];
			}
			else {
				for (unsigned int mu = 0; mu < 4; ++mu) even_inverse[site][mu] = odd_inverse[site][mu];
			}
		}
		eo->multiplyEvenEvenInverse(even_inverse);

		solution = even_inverse;
	} else {
		if (isOutputProcess()) std::cout << "PreconditionedBiCGStab::Not using even-odd preconditioning ... " << std::endl;
		solution = source;
		AlgebraUtils::gamma5(solution);
		reduced_dirac_vector_t tmp;
		dirac->multiply(tmp, solution);
		AlgebraUtils::gamma5(tmp);

		DiracOperator* squareDiracOperator = DiracOperator::getSquare(dirac);
		squareDiracOperator->setGamma5(true);

		res = biConjugateGradient->solve(squareDiracOperator, tmp, solution, initial_guess);
		delete squareDiracOperator;
	}

	lastSteps = biConjugateGradient->getLastSteps();
	return res;
}

} /* namespace Update */

