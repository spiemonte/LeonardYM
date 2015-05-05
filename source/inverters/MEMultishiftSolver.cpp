/*
 * MEMultishiftSolver.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#include "MEMultishiftSolver.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

MEMultishiftSolver::MEMultishiftSolver(real_t _epsilon, unsigned int _maxSteps) : MultishiftSolver(_epsilon, _maxSteps) { }

MEMultishiftSolver::~MEMultishiftSolver() { }

bool MEMultishiftSolver::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, std::vector<extended_dirac_vector_t>& original_solutions, const std::vector<real_t>& shifts) {
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

	long_real_t normResidual = AlgebraUtils::squaredNorm(residual);
	for (unsigned int i = 0; i < maxSteps; ++i) {
		dirac->multiplyAdd(tmp, p, p, *shift);
		long_real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));

#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*solution)[site][mu] = (*solution)[site][mu] + alpha*p[site][mu];
				residual[site][mu] = residual[site][mu] - alpha*tmp[site][mu];
			}
		}
		solution->updateHalo();
		residual.updateHalo();//TODO maybe not needed

		long_real_t error = AlgebraUtils::squaredNorm(residual);

		if (error < epsilon) {
#ifdef MULTISHIFTLOG
			std::cout << "Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
			break;
		}

		long_real_t beta = error / normResidual;

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
	std::vector<reduced_dirac_vector_t>::iterator previous_solution = solution;
	++solution;
	++shift;
	for (; shift != shifts.end(); ++solution, ++shift, ++previous_solution) {
		dirac->multiplyAdd(tmp, *previous_solution, *previous_solution, *shift);
		//First set the initial residual and p to source, result to the previous result
#pragma omp parallel for
		for (int site = 0; site < residual.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*solution)[site][mu] = (*previous_solution)[site][mu];
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
			long_real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));

#pragma omp parallel for
			for (int site = 0; site < residual.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					(*solution)[site][mu] = (*solution)[site][mu] + alpha*p[site][mu];
					residual[site][mu] = residual[site][mu] - alpha*tmp[site][mu];
				}
			}
			solution->updateHalo();
			residual.updateHalo();//TODO maybe not needed

			long_real_t error = AlgebraUtils::squaredNorm(residual);

			if (error < epsilon) {
#ifdef MULTISHIFTLOG
				if (isOutputProcess()) std::cout << "Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
				break;
			}
#ifdef FULLLOG
            else if (isOutputProcess()) std::cout << "Error " << error << " at step " << i << " for shift " << *shift << " and alpha: " << alpha << " and epsilon " << epsilon << std::endl;
#endif

			long_real_t beta = error / normResidual;

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

	if (!noproblem && isOutputProcess()) std::cout << "Failure in finding convergence in MultishiftSolver!" << std::endl;

	return noproblem;
}

} /* namespace Update */
