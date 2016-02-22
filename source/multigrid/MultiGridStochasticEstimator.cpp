#include "MultiGridStochasticEstimator.h"
#include "MultiGridProjector.h"

namespace Update {

MultiGridStochasticEstimator::MultiGridStochasticEstimator() : StochasticEstimator(), blockBasis(0), biMgSolver(0) { }

MultiGridStochasticEstimator::MultiGridStochasticEstimator(BlockBasis* _blockBasis) : StochasticEstimator(), blockBasis(_blockBasis), biMgSolver(0) { }

MultiGridStochasticEstimator::MultiGridStochasticEstimator(const MultiGridStochasticEstimator& toCopy) : StochasticEstimator(toCopy), blockBasis(toCopy.blockBasis), biMgSolver(new MultiGridBiConjugateGradientSolver()) { }

MultiGridStochasticEstimator::~MultiGridStochasticEstimator() {
	if (biMgSolver != 0) delete biMgSolver;
}

void MultiGridStochasticEstimator::getMultigridVectors(DiracOperator* dirac, extended_dirac_vector_t& mg_source, extended_dirac_vector_t& inverse_mg_source) {
	multigrid_vector_t mg_sr, mg_sol;

#pragma omp parallel for
	for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) {
#ifndef MULTITHREADING
		real_t realPart = (randomInteger() == 0 ? -1 : 1);
#endif
#ifdef MULTITHREADING
		real_t realPart = ((*randomInteger[omp_get_thread_num()])() == 0 ? -1 : 1);
#endif
		mg_sr[i] = std::complex<real_t>(realPart,0.);
	}

	MultiGridOperator* multiGridOperator = new MultiGridOperator();
	MultiGridProjector* multiGridProjector = new MultiGridProjector();
	multiGridOperator->setDiracOperator(dirac);

	for (int i = 0; i < blockBasis->size(); ++i) {
		multiGridOperator->addVector((*blockBasis)[i]);
		multiGridProjector->addVector((*blockBasis)[i]);
	}

	if (biMgSolver == 0) biMgSolver = new MultiGridBiConjugateGradientSolver();

	biMgSolver->setMaximumSteps(400);
	biMgSolver->setPrecision(0.0000000001);

	biMgSolver->solve(multiGridOperator, mg_sr, mg_sol);

	reduced_dirac_vector_t reduced_mg_source, reduced_inverse_mg_source;
	multiGridProjector->apply(reduced_mg_source,mg_sr);
	mg_source = reduced_mg_source;
	multiGridProjector->apply(reduced_inverse_mg_source,mg_sol);
	inverse_mg_source = reduced_inverse_mg_source;
	
	delete multiGridOperator;
	delete multiGridProjector;
}

void MultiGridStochasticEstimator::getResidualVectors(Solver* solver, DiracOperator* dirac, extended_dirac_vector_t& residual_source, extended_dirac_vector_t& inverse_residual_source) {
	multigrid_vector_t mg_sr, mg_sol;

	this->generateRandomNoise(residual_source);
	reduced_dirac_vector_t reduced_residual_source = residual_source, mg_inverse, full_inverse;

	MultiGridOperator* multiGridOperator = new MultiGridOperator();
	MultiGridProjector* multiGridProjector = new MultiGridProjector();
	multiGridOperator->setDiracOperator(dirac);

	for (int i = 0; i < blockBasis->size(); ++i) {
		multiGridOperator->addVector((*blockBasis)[i]);
		multiGridProjector->addVector((*blockBasis)[i]);
	}

	if (biMgSolver == 0) biMgSolver = new MultiGridBiConjugateGradientSolver();

	biMgSolver->setMaximumSteps(400);
	biMgSolver->setPrecision(0.0000000001);

	multiGridProjector->apply(mg_sr, reduced_residual_source);

	biMgSolver->solve(multiGridOperator, mg_sr, mg_sol);

	multiGridProjector->apply(mg_inverse,mg_sol);
	solver->solve(dirac, reduced_residual_source, full_inverse);

#pragma omp parallel for
	for (int site = 0; site < full_inverse.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			full_inverse[site][mu] = full_inverse[site][mu] - mg_inverse[site][mu];
		}
	}

	inverse_residual_source = full_inverse;

	delete multiGridOperator;
	delete multiGridProjector;
}

void MultiGridStochasticEstimator::setBlockBasis(BlockBasis* _blockBasis) {
	blockBasis = _blockBasis;
}

}

