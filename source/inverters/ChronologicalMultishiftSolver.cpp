/*
 * ChronologicalMultishiftSolver.cpp
 *
 *  Created on: Jul 24, 2012
 *      Author: spiem_01
 */

#include "ChronologicalMultishiftSolver.h"
#include "AlgebraUtils.h"
#define FULLLOG

namespace Update {

ChronologicalMultishiftSolver::ChronologicalMultishiftSolver(double _epsilon, unsigned int _maxSteps) : MultishiftSolver(_epsilon, _maxSteps) { }

ChronologicalMultishiftSolver::~ChronologicalMultishiftSolver() { }

bool ChronologicalMultishiftSolver::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, std::vector<extended_dirac_vector_t>& original_solutions, const std::vector<real_t>& shifts) {
	//We work with reduced halos
	reduced_dirac_vector_t source = original_source;
	std::vector<reduced_dirac_vector_t> solutions(original_solutions.size());

	
	//First solve the linear system of the first shift (supposed to be the easiest, i.e. s[0] > s[i])
	std::vector<reduced_dirac_vector_t>::iterator solution = solutions.begin();
	std::vector<real_t>::const_iterator shift = shifts.begin();
	//The flag for the errors
	bool noproblem = true;
	//First set the initial residual and p to source, result to zero
#pragma omp parallel for
	for (int site = 0; site < residual.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero((*solution)[site][mu]);
			p[site][mu] = source[site][mu];
			residual[site][mu] = source[site][mu];
		}
	}
	solution->updateHalo();
	p.updateHalo();
	residual.updateHalo();

	real_t normResidual = AlgebraUtils::squaredNorm(residual);
	for (unsigned int i = 0; i < maxSteps; ++i) {
		dirac->multiplyAdd(tmp, p, p, *shift);
		real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));


#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*solution)[site][mu] = (*solution)[site][mu] + alpha*p[site][mu];
				residual[site][mu] = residual[site][mu] - alpha*tmp[site][mu];
			}
		}
		solution->updateHalo();
		residual.updateHalo();

		real_t error = AlgebraUtils::squaredNorm(residual);

		if (error < epsilon) {
#ifdef FULLLOG
			if (isOutputProcess()) std::cout << "Chronological Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
			break;
		}

		real_t beta = error / normResidual;


#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = residual[site][mu] + beta*p[site][mu];
			}
		}
		p.updateHalo();

		normResidual = error;

		if (i == maxSteps - 1) noproblem = false;
	}


	//Now we have the first solution, we can use the previous solution as initial guess for the others
	//We use also the old chronological solution as initial guess, minimizing the norm
	std::vector<reduced_dirac_vector_t>::iterator previous_solution = solution;
	++solution;
	++shift;
	for (; shift != shifts.end(); ++solution, ++shift, ++previous_solution) {
		//We use now p and tmp as temporary storage
		dirac->multiplyAdd(p, *previous_solution, *previous_solution, *shift);
		dirac->multiplyAdd(tmp, *solution, *solution, *shift);
		//Now we solve the CG for p and tmp
		long_real_t pApA = AlgebraUtils::squaredNorm(p);
		long_real_t hAhA = AlgebraUtils::squaredNorm(tmp);
		std::complex<long_real_t> pAb = AlgebraUtils::dot(p,source);
		std::complex<long_real_t> hAb = AlgebraUtils::dot(tmp,source);
		std::complex<long_real_t> hApA = AlgebraUtils::dot(tmp,p);
		long_real_t c1R = real((-(hAhA*pAb) + hApA*conj(hAb) + hAb*conj(hApA) - hAhA*conj(pAb))/(std::complex<long_real_t>(2.,0.)*(-(hAhA*pApA) + hApA*conj(hApA))));
		long_real_t c1I = real((std::complex<long_real_t>(0,-0.5)*(-(hAhA*pAb) - hApA*conj(hAb) + hAb*conj(hApA) + hAhA*conj(pAb)))/(-(hAhA*pApA) + hApA*conj(hApA)));
		long_real_t c2R = real((hApA*pAb - hAb*pApA - pApA*conj(hAb) + conj(hApA)*conj(pAb))/(std::complex<long_real_t>(2.,0.)*(-(hAhA*pApA) + hApA*conj(hApA))));
		long_real_t c2I = real((std::complex<long_real_t>(0,-0.5)*(hApA*pAb - hAb*pApA + pApA*conj(hAb) - conj(hApA)*conj(pAb)))/(-(hAhA*pApA) + hApA*conj(hApA)));
		std::complex<real_t> c1(c1R,c1I);
		std::complex<real_t> c2(c2R,c2I);

		//We are ready for the new guess
#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				initial_guess[site][mu] = c1*(*previous_solution)[site][mu]+c2*(*solution)[site][mu];
			}
		}
		initial_guess.updateHalo();

		/*dirac->multiplyAdd(tmp, initial_guess, initial_guess, *shift);
		long_real_t test = 0.;
		for (unsigned int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				test += real((tmp[site][mu] - source[site][mu]).dot(tmp[site][mu] - source[site][mu]));
			}
		}
		std::cout << "Dopo: " << test << std::endl;
		dirac->multiplyAdd(tmp, *previous_solution, *previous_solution, *shift);
		test = 0.;
		for (unsigned int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				test += real((tmp[site][mu] - source[site][mu]).dot(tmp[site][mu] - source[site][mu]));
			}
		}
		std::cout << "Prima: " << test << std::endl;*/

		//Now we start as usual with the biCGStab
		dirac->multiplyAdd(tmp, initial_guess, initial_guess, *shift);

		//std::cout << "Vediamo un po': " << AlgebraUtils::squaredNorm(initial_guess) << " " << AlgebraUtils::squaredNorm(*previous_solution) << " " << c1 << " " << c2 << std::endl;


		//First set the initial residual and p to source, result to the previous result
#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*solution)[site][mu] = initial_guess[site][mu];
				residual[site][mu] = source[site][mu] - tmp[site][mu];
				p[site][mu] = source[site][mu] - tmp[site][mu];
			}
		}
		solution->updateHalo();
		residual.updateHalo();
		p.updateHalo();

		normResidual = AlgebraUtils::squaredNorm(residual);
		for (unsigned int i = 0; i < maxSteps; ++i) {
			dirac->multiplyAdd(tmp, p, p, *shift);
			real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));

#pragma omp parallel for
			for (int site = 0; site < residual.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					(*solution)[site][mu] = (*solution)[site][mu] + alpha*p[site][mu];
					residual[site][mu] = residual[site][mu] - alpha*tmp[site][mu];
				}
			}
			solution->updateHalo();
			residual.updateHalo();//TODO maybe not needed

			real_t error = AlgebraUtils::squaredNorm(residual);

			if (error < epsilon) {
#ifdef FULLLOG
				if (isOutputProcess()) std::cout << "Chronological Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
				break;
			}

			real_t beta = error / normResidual;

#pragma omp parallel for
			for (int site = 0; site < residual.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					p[site][mu] = residual[site][mu] + beta*p[site][mu];
				}
			}
			p.updateHalo();

			normResidual = error;

			if (i == maxSteps - 1) noproblem = false;
		}

	}

	for (unsigned int index = 0; index < shifts.size(); ++index) {
		original_solutions[index] = solutions[index];
	}

	if (!noproblem) std::cout << "Failure in finding convergence in MultishiftSolver!" << std::endl;

	return noproblem;
}

} /* namespace Update */
