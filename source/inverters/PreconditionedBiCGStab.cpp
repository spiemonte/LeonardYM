#include "PreconditionedBiCGStab.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_operators/SquareEvenOddImprovedDiracWilsonOperator.h"

namespace Update {

PreconditionedBiCGStab::PreconditionedBiCGStab(const PreconditionedBiCGStab& toCopy) : Solver(toCopy), useEvenOddPreconditioning(true), biConjugateGradient(0) { }

PreconditionedBiCGStab::PreconditionedBiCGStab() : Solver(), useEvenOddPreconditioning(true), biConjugateGradient(0) { }

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

bool PreconditionedBiCGStab::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, DiracOperator* , reduced_dirac_vector_t const* initial_guess) {
	typedef reduced_dirac_vector_t::Layout Layout;
	bool res = false;

	if (biConjugateGradient == 0) {
		biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setMaximumSteps(this->getMaximumSteps());
		biConjugateGradient->setPrecision(this->getPrecision());
	}

	if (dynamic_cast<ImprovedDiracWilsonOperator*>(dirac) && (Layout::glob_x % 2 == 0) && (Layout::glob_y % 2 == 0) && (Layout::glob_z % 2 == 0) && (Layout::glob_t % 2 == 0) && useEvenOddPreconditioning) {
		ImprovedDiracWilsonOperator* dirac_improved = dynamic_cast<ImprovedDiracWilsonOperator*>(dirac);
		
		SquareEvenOddImprovedDiracWilsonOperator* eo_sq = new SquareEvenOddImprovedDiracWilsonOperator();
		eo_sq->setKappa(dirac_improved->getKappa());
		eo_sq->setCSW(dirac_improved->getCSW());
		eo_sq->setLattice(*dirac->getLattice());

		EvenOddImprovedDiracWilsonOperator* eo = new EvenOddImprovedDiracWilsonOperator();
		eo->setKappa(dirac_improved->getKappa());
		eo->setCSW(dirac_improved->getCSW());
		eo->setLattice(*dirac->getLattice());
	
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

		delete eo;
		delete eo_sq;
	} else {
		if (isOutputProcess()) std::cout << "PreconditionedBiCGStab::Not using even-odd preconditioning ... " << std::endl;

		res = biConjugateGradient->solve(dirac, source, solution, initial_guess);
	}

	lastSteps = biConjugateGradient->getLastSteps();
	return res;
}

void PreconditionedBiCGStab::setUseEvenOddPreconditioning(bool _useEvenOddPreconditioning) {
	useEvenOddPreconditioning = _useEvenOddPreconditioning;	
}

} /* namespace Update */

