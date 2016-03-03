/*
 * MEMultishiftSolver.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#include "MultiGridMEMultishiftSolver.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "GMRESR.h"

namespace Update {

MultiGridMEMultishiftSolver::MultiGridMEMultishiftSolver(int basisDimension, const std::vector<unsigned int>& _blockSize, BlockDiracOperator* _blackBlockDiracOperator, BlockDiracOperator* _redBlockDiracOperator, real_t _epsilon, unsigned int _maxSteps) : MultishiftSolver(_epsilon, _maxSteps), MultiGridSolver(basisDimension, _blockSize, _blackBlockDiracOperator, _redBlockDiracOperator) { }

MultiGridMEMultishiftSolver::MultiGridMEMultishiftSolver(const MultiGridMEMultishiftSolver& toCopy) : MultishiftSolver(toCopy), MultiGridSolver(toCopy) { }

MultiGridMEMultishiftSolver::~MultiGridMEMultishiftSolver() { }

bool MultiGridMEMultishiftSolver::solve(DiracOperator* diracSquare, const extended_dirac_vector_t& original_source, std::vector<extended_dirac_vector_t>& original_solutions, const std::vector<real_t>& shifts) {
	//We work with reduced halos
	reduced_dirac_vector_t source = original_source;
	std::vector<reduced_dirac_vector_t> solutions(original_solutions.size());
	std::vector<reduced_dirac_vector_t> square_root_solutions(original_solutions.size());

	DiracOperator* dirac = DiracOperator::getSquareRoot(diracSquare);
	dirac->setGamma5(false);
	TwistedDiracOperator* twistedDirac = new TwistedDiracOperator();
	twistedDirac->setDiracOperator(dirac);
	
	GMRESR* gmresr = new GMRESR();
	
	//First solve the linear system of the first shift (supposed to be the easiest, i.e. s[0] > s[i])
	std::vector<reduced_dirac_vector_t>::iterator solution = solutions.begin();
	std::vector<reduced_dirac_vector_t>::iterator square_root_solution = square_root_solutions.begin();
	std::vector<real_t>::const_iterator shift = shifts.begin();

	//consider first the negative and large positive shifts
	while (shift != shifts.end() && (*shift < 0 || *shift > 0.001)) {
		if (isOutputProcess()) std::cout << "MultiGridMEMultishiftSolver::Inverting shift " << *shift << std::endl;
		//Solve it separately
		///
		reduced_dirac_vector_t tmp;
		
#pragma omp parallel for
		for (int site = 0; site < r.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_zero((*solution)[site][mu]);
				p[site][mu] = source[site][mu];
				r[site][mu] = source[site][mu];
			}
		}

		long_real_t normResidual = AlgebraUtils::squaredNorm(r);
		for (unsigned int i = 0; i < MultishiftSolver::maxSteps; ++i) {
			diracSquare->multiplyAdd(tmp, p, p, *shift);
			long_real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));

#pragma omp parallel for
			for (int site = 0; site < r.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					(*solution)[site][mu] = (*solution)[site][mu] + alpha*p[site][mu];
					r[site][mu] = r[site][mu] - alpha*tmp[site][mu];
				}
			}

			long_real_t error = AlgebraUtils::squaredNorm(r);

			if (error < epsilon) {
#ifdef MULTISHIFTLOG
				std::cout << "Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
				break;
			}

			long_real_t beta = error / normResidual;

#pragma omp parallel for
			for (int site = 0; site < r.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					p[site][mu] = r[site][mu] + beta*p[site][mu];
				}
			}

			normResidual = error;
		}

		///
		++shift;
		++solution;
		++square_root_solution;
	}

	if (shift == shifts.end()) {
		for (unsigned int index = 0; index < shifts.size(); ++index) {
			original_solutions[index] = solutions[index];
		}
		return true;
	}

	AlgebraUtils::gamma5(source);
	
	//We solve the first linear equation
	if (isOutputProcess()) std::cout << "MultiGridMEMultishiftSolver::Inverting shift " << *shift << std::endl;
	twistedDirac->setTwist(-sqrt(*shift));
	redBlockDiracOperator->setTwist(-sqrt(*shift));
	blackBlockDiracOperator->setTwist(-sqrt(*shift));
	this->solve(twistedDirac, source, *square_root_solution);
	
	AlgebraUtils::gamma5(*square_root_solution);
	
	twistedDirac->setTwist(sqrt(*shift));
	redBlockDiracOperator->setTwist(sqrt(*shift));
	blackBlockDiracOperator->setTwist(sqrt(*shift));
	this->solve(twistedDirac, *square_root_solution, *solution);
	
	//restore
	AlgebraUtils::gamma5(*square_root_solution);
	//AlgebraUtils::gamma5(source);
	//AlgebraUtils::gamma5(*solution);

	//reduced_dirac_vector_t tmp;
	//diracSquare->multiplyAdd(tmp, *solution, *solution, *shift);
	//twistedDirac->multiply(tmp, source);
	//std::cout << "Vediamo di che morte dobbiamo morire: " << AlgebraUtils::differenceNorm(tmp, source) << " " << (*solution)[5][2][3] << " " << AlgebraUtils::squaredNorm(*solution) << " " << sqrt(*shift) << " " << *shift << std::endl;


	//Now we have the first solution, we can use the previous solution as initial guess for the others
	std::vector<reduced_dirac_vector_t>::iterator previous_solution = solution;
	std::vector<reduced_dirac_vector_t>::iterator previous_square_root_solution = square_root_solution;
	++solution;
	++square_root_solution;
	++shift;
	for (; shift != shifts.end(); ++solution, ++shift, ++previous_solution, ++previous_square_root_solution) {
		if (isOutputProcess()) std::cout << "MultiGridMEMultishiftSolver::Inverting shift " << *shift << std::endl;
		twistedDirac->setTwist(-sqrt(*shift));
		redBlockDiracOperator->setTwist(-sqrt(*shift));
		blackBlockDiracOperator->setTwist(-sqrt(*shift));
		this->solve(twistedDirac, source, *square_root_solution, &(*previous_square_root_solution));

		AlgebraUtils::gamma5(*square_root_solution);
		
		twistedDirac->setTwist(sqrt(*shift));
		redBlockDiracOperator->setTwist(sqrt(*shift));
		blackBlockDiracOperator->setTwist(sqrt(*shift));
		this->solve(twistedDirac, *square_root_solution, *solution, &(*previous_solution));
		
		//restore
		AlgebraUtils::gamma5(*square_root_solution);
	}

	for (unsigned int index = 0; index < shifts.size(); ++index) {
		original_solutions[index] = solutions[index];
	}

	delete  twistedDirac;
	delete gmresr;
	delete dirac;
	redBlockDiracOperator->setTwist(0.);
	blackBlockDiracOperator->setTwist(0.);

	return true;
}

} /* namespace Update */
