#include "MultiGridSolver.h"
#include "MultiGridProjector.h"
#include "algebra_utils/AlgebraUtils.h"
#include "inverters/GMRESR.h"

namespace Update {

MultiGridSolver::MultiGridSolver(int basisDimension, const std::vector<unsigned int>& _blockSize, BlockDiracOperator* _blackBlockDiracOperator, BlockDiracOperator* _redBlockDiracOperator) : Solver("MultiGridSolver"), blockBasis(basisDimension), blockSize(_blockSize), blackBlockDiracOperator(_blackBlockDiracOperator), redBlockDiracOperator(_redBlockDiracOperator), K(0), biMgSolver(new MultiGridBiConjugateGradientSolver()), SAPIterantions(7), SAPMaxSteps(100), SAPPrecision(0.00001), GMRESIterations(300), GMRESPrecision(0.0000000001), BiMGIterations(35), BiMGPrecision(0.00000000001) { }

bool MultiGridSolver::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution, reduced_dirac_vector_t const* initial_guess) {
	blackBlockDiracOperator->setLattice(*dirac->getLattice());
	blackBlockDiracOperator->setGamma5(false);
	blackBlockDiracOperator->setBlockSize(blockSize);
		
	redBlockDiracOperator->setLattice(*dirac->getLattice());
	redBlockDiracOperator->setGamma5(false);
	redBlockDiracOperator->setBlockSize(blockSize);
		
	K = new ComplementBlockDiracOperator(dirac, redBlockDiracOperator, blackBlockDiracOperator);
	K->setMaximumSteps(SAPMaxSteps);
	K->setBlockSize(blockSize);

	SAPPreconditioner *preconditioner = new SAPPreconditioner(dirac,K);
	preconditioner->setSteps(SAPIterantions);
	preconditioner->setPrecision(SAPPrecision);

	MultiGridOperator* multiGridOperator = new MultiGridOperator();
	MultiGridProjector* multiGridProjector = new MultiGridProjector();
	multiGridOperator->setDiracOperator(dirac);

	for (int i = 0; i < blockBasis.size(); ++i) {
		multiGridOperator->addVector(blockBasis[i]);
		multiGridProjector->addVector(blockBasis[i]);
	}

	biMgSolver->setMaximumSteps(BiMGIterations);
	biMgSolver->setPrecision(BiMGPrecision);

	multigrid_vector_t solution_hat, source_hat;


	int conjugateSpaceDimension = 15;

	if (initial_guess == 0) {
		long_real_t normSource = AlgebraUtils::squaredNorm(source);
		if (normSource > precision) {
			solution = source;
		}
		else {
			AlgebraUtils::setToZero(solution);
		}
	}
	else {
		solution = *initial_guess;
	}

	dirac->multiply(r,solution);
#pragma omp parallel for
	for (int site = 0; site < solution.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) r[site][mu] = source[site][mu] - r[site][mu];
	}
	
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);

	for (unsigned int k = 0; k < maxSteps; ++k) {
		{
			multiGridProjector->apply(source_hat,r);
			biMgSolver->solve(multiGridOperator, source_hat, solution_hat);

			multiGridProjector->apply(mg_inverse,solution_hat);

			dirac->multiply(source_sap,mg_inverse);
#pragma omp parallel for
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					source_sap[site][mu] = r[site][mu] - source_sap[site][mu];
				}
			}
			
			preconditioner->multiply(u[(k % conjugateSpaceDimension)],source_sap);
#pragma omp parallel for
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu] + mg_inverse[site][mu];
				}
			}

		}
		//u[(k % conjugateSpaceDimension)] = r;

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
		lastError = error;
		if (error < precision) {
			maxSteps = k;
			break;
		}
		else if (isOutputProcess()) {
			std::cout << "MultiGridSolver::Residual norm at step " << k << ": " << error << std::endl;
		}

	}

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	if (isOutputProcess()) std::cout << "MultiGridSolver::Multigrid inversion done in: " << (elapsed) << " s."<< std::endl;

	delete K;
	delete preconditioner;
	delete multiGridOperator;
	delete multiGridProjector;

	return true;
	
}

void MultiGridSolver::initializeBasis(DiracOperator* dirac) {
	if (MultiGridVectorLayout::xBlockSize != blockSize[0] || MultiGridVectorLayout::yBlockSize != blockSize[1] || MultiGridVectorLayout::zBlockSize != blockSize[2] || MultiGridVectorLayout::tBlockSize != blockSize[3]) {
		MultiGridVectorLayout::xBlockSize = blockSize[0];
		MultiGridVectorLayout::yBlockSize = blockSize[1];
		MultiGridVectorLayout::zBlockSize = blockSize[2];
		MultiGridVectorLayout::tBlockSize = blockSize[3];
		MultiGridVectorLayout::initialize();
	}

	reduced_dirac_vector_t zeroVector, tmp;
	AlgebraUtils::setToZero(zeroVector);
	reduced_dirac_vector_t randomVector;
	GMRESR* gmres_inverter = new GMRESR();
	gmres_inverter->setPrecision(GMRESPrecision);
	gmres_inverter->setMaximumSteps(GMRESIterations);

	blackBlockDiracOperator->setLattice(*dirac->getLattice());
	blackBlockDiracOperator->setGamma5(false);
	blackBlockDiracOperator->setBlockSize(blockSize);
		
	redBlockDiracOperator->setLattice(*dirac->getLattice());
	redBlockDiracOperator->setGamma5(false);
	redBlockDiracOperator->setBlockSize(blockSize);
		
	K = new ComplementBlockDiracOperator(dirac, redBlockDiracOperator, blackBlockDiracOperator);
	K->setMaximumSteps(205);
	K->setBlockSize(blockSize);

	SAPPreconditioner *preconditioner = new SAPPreconditioner(dirac,K);
	preconditioner->setSteps(SAPIterantions);
	preconditioner->setPrecision(SAPPrecision);

	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_REALTIME, &start);
		
	for (int i = 0; i < blockBasis.size(); i += 1) {
		//We start with a random vector
		AlgebraUtils::generateRandomVector(randomVector);
		AlgebraUtils::normalize(randomVector);

		//We give random vector as initial guess to solve the omogeneous system
		gmres_inverter->solve(dirac, zeroVector, tmp, preconditioner, &randomVector);
		gmres_inverter->solve(dirac, tmp, blockBasis[i], preconditioner);
		//gmres_inverter->solve(dirac, blockBasis[i], tmp, preconditioner);
		//gmres_inverter->solve(dirac, tmp, blockBasis[i], preconditioner);
		
		dirac->multiply(randomVector, blockBasis[i]);
		
		long_real_t deficit = AlgebraUtils::squaredNorm(randomVector)/AlgebraUtils::squaredNorm(blockBasis[i]);
		if (isOutputProcess()) std::cout << "MultiGridSolver::Deficit for the basis vector " << i << " " << sqrt(deficit) << std::endl;
	}

	blockBasis.orthogonalize();

	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	if (isOutputProcess()) std::cout << "MultiGridSolver::Multigrid basis constructed in: " << (elapsed) << " s."<< std::endl;

	delete gmres_inverter;
	delete K;
	delete preconditioner;
}

void MultiGridSolver::setSAPIterations(int _SAPIterantions) {
	SAPIterantions = _SAPIterantions;
}

int MultiGridSolver::getSAPIterations() const {
	return SAPIterantions;
}

void MultiGridSolver::setSAPMaxSteps(int _SAPMaxSteps) {
	SAPMaxSteps = _SAPMaxSteps;
}

int MultiGridSolver::getSAPMaxSteps() const {
	return SAPMaxSteps;
}

void MultiGridSolver::setSAPPrecision(real_t _SAPPrecision) {
	SAPPrecision = _SAPPrecision;
}

real_t MultiGridSolver::getSAPPrecision() const {
	return SAPPrecision;
}

void MultiGridSolver::setGMRESIterations(int _GMRESIterations) {
	GMRESIterations = _GMRESIterations;
}

int MultiGridSolver::getGMRESIterations() const {
	return GMRESIterations;
}

void MultiGridSolver::setGMRESPrecision(real_t _GMRESPrecision) {
	GMRESPrecision = _GMRESPrecision;
}

real_t MultiGridSolver::getGMRESPrecision() const {
	return GMRESPrecision;
}

void MultiGridSolver::setBasisDimension(unsigned int dim) {
	blockBasis.setBasisDimension(dim);
}

unsigned int MultiGridSolver::getBasisDimension() const {
	return blockBasis.size();
}
	
BlockBasis* MultiGridSolver::getBasis() {
	return &blockBasis;
}

}

