#include "ExactOverlapOperator.h"
#include "hmc_forces/OverlapFermionForce.h"
#include "fermion_measurements/DiracEigenSolver.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

ExactOverlapOperator::ExactOverlapOperator() : OverlapOperator(), diracEigenSolver(new DiracEigenSolver()), recomputeEigenvalues(true), numberOfEigenvalues(100) { }

ExactOverlapOperator::ExactOverlapOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : OverlapOperator(_lattice, _kappa, _gamma5), diracEigenSolver(0), recomputeEigenvalues(true), numberOfEigenvalues(100) { }

ExactOverlapOperator::ExactOverlapOperator(const ExactOverlapOperator& copy) : OverlapOperator(copy.lattice, copy.kappa, copy.gamma5), diracEigenSolver(new DiracEigenSolver(*copy.diracEigenSolver)), recomputeEigenvalues(copy.recomputeEigenvalues), computed_eigenvalues(copy.computed_eigenvalues), computed_eigenvectors(copy.computed_eigenvectors), numberOfEigenvalues(copy.numberOfEigenvalues) {}

ExactOverlapOperator::~ExactOverlapOperator() {
	if (diracEigenSolver) delete diracEigenSolver;
}

void ExactOverlapOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	if (recomputeEigenvalues) this->computeEigenvalues();
	squareDiracWilsonOperator.setGamma5(true);
	diracWilsonOperator.setGamma5(true);

	AlgebraUtils::setToZero(output);

	tmp1 = input;

	for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
		std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(computed_eigenvectors[i], input));
		std::complex<real_t> sign_proj = proj;
		if (computed_eigenvalues[i] < 0.) sign_proj = -proj;

#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += sign_proj*computed_eigenvectors[i][site][mu];
				tmp1[site][mu] -= proj*computed_eigenvectors[i][site][mu];
			}
		}
	}

	squareRootApproximation.evaluate(&squareDiracWilsonOperator, tmp2, tmp1);
	diracWilsonOperator.multiply(tmp1, tmp2);
	
	if (gamma5) {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*input[site][0] + (0.5-mass/2.)*(tmp1[site][0] + output[site][0]);
			output[site][1] = (0.5+mass/2.)*input[site][1] + (0.5-mass/2.)*(tmp1[site][1] + output[site][1]);
			output[site][2] = -(0.5+mass/2.)*input[site][2] + (0.5-mass/2.)*(tmp1[site][2] + output[site][2]);
			output[site][3] = -(0.5+mass/2.)*input[site][3] + (0.5-mass/2.)*(tmp1[site][3] + output[site][3]);
		}
	}
	else {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*input[site][0] + (0.5-mass/2.)*(tmp1[site][0] + output[site][0]);
			output[site][1] = (0.5+mass/2.)*input[site][1] + (0.5-mass/2.)*(tmp1[site][1] + output[site][1]);
			output[site][2] = (0.5+mass/2.)*input[site][2] - (0.5-mass/2.)*(tmp1[site][2] + output[site][2]);
			output[site][3] = (0.5+mass/2.)*input[site][3] - (0.5-mass/2.)*(tmp1[site][3] + output[site][3]);
		}
	}

	output.updateHalo();
}

void ExactOverlapOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	if (recomputeEigenvalues) this->computeEigenvalues();
	squareDiracWilsonOperator.setGamma5(true);
	diracWilsonOperator.setGamma5(true);

	AlgebraUtils::setToZero(output);

	tmp1 = vector1;

	for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
		std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(computed_eigenvectors[i], vector1));
		std::complex<real_t> sign_proj = proj;
		if (computed_eigenvalues[i] < 0.) sign_proj = -proj;

#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += sign_proj*computed_eigenvectors[i][site][mu];
				tmp1[site][mu] -= proj*computed_eigenvectors[i][site][mu];
			}
		}
	}
	
	squareRootApproximation.evaluate(&squareDiracWilsonOperator, tmp2, tmp1);
	diracWilsonOperator.multiply(tmp1, tmp2);
	
	if (gamma5) {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*vector1[site][0] + (0.5-mass/2.)*(tmp1[site][0] + output[site][0]) + alpha*vector2[site][0];
			output[site][1] = (0.5+mass/2.)*vector1[site][1] + (0.5-mass/2.)*(tmp1[site][1] + output[site][1]) + alpha*vector2[site][1];
			output[site][2] = -(0.5+mass/2.)*vector1[site][2] + (0.5-mass/2.)*(tmp1[site][2] + output[site][2]) + alpha*vector2[site][2];
			output[site][3] = -(0.5+mass/2.)*vector1[site][3] + (0.5-mass/2.)*(tmp1[site][3] + output[site][3]) + alpha*vector2[site][3];
		}
	}
	else {
#pragma omp parallel for
		for (int site = 0; site < output.localsize; ++site) {
			output[site][0] = (0.5+mass/2.)*vector1[site][0] + (0.5-mass/2.)*(tmp1[site][0] + output[site][0]) + alpha*vector2[site][0];
			output[site][1] = (0.5+mass/2.)*vector1[site][1] + (0.5-mass/2.)*(tmp1[site][1] + output[site][1]) + alpha*vector2[site][1];
			output[site][2] = (0.5+mass/2.)*vector1[site][2] - (0.5-mass/2.)*(tmp1[site][2] + output[site][2]) + alpha*vector2[site][2];
			output[site][3] = (0.5+mass/2.)*vector1[site][3] - (0.5-mass/2.)*(tmp1[site][3] + output[site][3]) + alpha*vector2[site][3];
		}
	}

	output.updateHalo();
}

void ExactOverlapOperator::computeEigenvalues() {
	if (diracEigenSolver == 0 || !this->checkEigenvalues()) {
		std::vector< std::complex<real_t> > squared_eigenvalues;

		diracEigenSolver->maximumEigenvalues(&squareDiracWilsonOperator, squared_eigenvalues, computed_eigenvectors, numberOfEigenvalues, SmallestReal);
	
		computed_eigenvalues.resize(squared_eigenvalues.size());

		for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
			diracWilsonOperator.multiply(tmp1, computed_eigenvectors[i]);
			computed_eigenvalues[i] = static_cast< real_t >(real(AlgebraUtils::dot(computed_eigenvectors[i],tmp1)));
			if (computed_eigenvalues[i] < 0.) computed_eigenvalues[i] = -sqrt(real(squared_eigenvalues[i]));
			else computed_eigenvalues[i] = sqrt(real(squared_eigenvalues[i]));

			long_real_t convergence = 0.;
#pragma omp parallel for reduction(+:convergence)
			for (int site = 0; site < computed_eigenvectors[i].localsize; ++site) {
				real_t sum = 0.;
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (unsigned int c = 0; c < diracVectorLength; ++c) {
					sum += std::abs(tmp1[site][mu][c]-computed_eigenvalues[i]*computed_eigenvectors[i][site][mu][c]);
					}
				}
				convergence += sum;
			}
			reduceAllSum(convergence);

			if (isOutputProcess()) std::cout << "ExactOverlapOperator::precision of eigenvalue " << i << " is " << convergence << ", lambda = " << computed_eigenvalues[i] << std::endl;
		}

		recomputeEigenvalues = false;
	}
}

bool ExactOverlapOperator::checkEigenvalues() {
	if (computed_eigenvectors.size() == 0) return false;
 
	for (unsigned int i = 0; i < computed_eigenvectors.size(); ++i) {
		diracWilsonOperator.multiply(tmp1, computed_eigenvectors[i]);
		
		long_real_t convergence = 0.;
#pragma omp parallel for reduction(+:convergence)
		for (int site = 0; site < computed_eigenvectors[i].localsize; ++site) {
			real_t sum = 0.;
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int c = 0; c < diracVectorLength; ++c) {
					sum += std::abs(tmp1[site][mu][c]-computed_eigenvalues[i]*computed_eigenvectors[i][site][mu][c]);
				}
			}
			convergence += sum;
		}
		reduceAllSum(convergence);
	
		std::cout << "ExactOverlapOperator::precision of eigenvalue " << i << " is " << convergence << ", lambda = " << computed_eigenvalues[i] << std::endl;

		if (convergence > diracEigenSolver->getTolerance()) {
			return false;
		}
	}

	return true;
}

void ExactOverlapOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracWilsonOperator.setKappa(_kappa);
	squareDiracWilsonOperator.setKappa(_kappa);
	recomputeEigenvalues = true;
}

void ExactOverlapOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	diracWilsonOperator.setLattice(_lattice);
	squareDiracWilsonOperator.setLattice(_lattice);
	recomputeEigenvalues = true;
}

DiracEigenSolver* ExactOverlapOperator::getDiracEigenSolver() const {
	return diracEigenSolver;
}

void ExactOverlapOperator::setNumberOfEigenvalues(unsigned int _numberOfEigenvalues) {
	numberOfEigenvalues = _numberOfEigenvalues;
}

unsigned int ExactOverlapOperator::getNumberOfEigenvalues() const {
	return numberOfEigenvalues;
}

} /* namespace Update */
