/*
 * AugConjugateGradientSolver.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: spiem_01
 */

#include "AugConjugateGradientSolver.h"

namespace Update {

AugConjugateGradient::AugConjugateGradient() { }

AugConjugateGradient::~AugConjugateGradient() { }

bool AugConjugateGradient::solve(DiracOperator* dirac, const reduced_dirac_vector_t& source, reduced_dirac_vector_t& solution) {
	solution = source;
	dirac->multiply(tmp,solution);

#pragma omp parallel for
	for (int site = 0; site < residual.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			residual[site][mu] = source[site][mu] - tmp[site][mu];
		}
	}

	for (int j = 0; j < m; ++j) {
		std::complex<real_t> rproj = AlgebraUtils::dot(residual, w[j]);
#pragma omp parallel for
		for (int site = 0; site < residual.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + (rproj/wproj[j])*w[j][site][mu];
				residual[site][mu] = residual[site][mu] - (rproj/wproj[j])*dw[j][site][mu];
			}
		}
	}

	zeta = residual;

	for (int j = 0; j < m; ++j) {
		std::complex<real_t> rproj = AlgebraUtils::dot(zeta, dw[j]);
#pragma omp parallel for
		for (int site = 0; site < residual.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				zeta[site][mu] = zeta[site][mu] - (rproj/wproj[j])*w[j][site][mu];
			}
		}
	}

	p = zeta;

	for (int k = 0; k < maxSteps; ++k) {
		dirac->multiply(tmp,p);

		std::complex<real_t> rz = AlgebraUtils::dot(residual, zeta);
		std::complex<real_t> pdp = AlgebraUtils::dot(p, tmp);
		std::complex<real_t> alpha = rz/pdp;

#pragma omp parallel for
		for (int site = 0; site < residual.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = solution[site][mu] + alpha*p[site][mu];
				residual[site][mu] = residual[site][mu] - alpha*tmp[site][mu];
			}
		}

		real_t residualNorm = AlgebraUtils::squaredNorm(residual);
		if (residualNorm < epsilon) {
			if (isOutputProcess()) std::cout << "AugConjugateGradient::Convergence in " << j << " steps" << std::endl;
			return true;
		}

		std::complex<real_t> mu = AlgebraUtils::dot(residual,dw[m-1])/wproj[m-1];

#pragma omp parallel for
		for (int site = 0; site < residual.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				zeta[site][mu] = residual[site][mu] - mu*dw[m-1][site][mu];
			}
		}

		std::complex<real_t> beta = AlgebraUtils::dot(residual,zeta)/rz;

#pragma omp parallel for
		for (int site = 0; site < residual.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				p[site][mu] = zeta[site][mu] + beta*p[site][mu];
			}
		}
	}

	if (isOutputProcess()) std::cout << "AugConjugateGradient::Failure in finding convergence in " << maxSteps << " steps" << std::endl;
	return false;
}

void AugConjugateGradient::setPrecision(double _epsilon) {
	epsilon = _epsilon;
}

double AugConjugateGradient::getPrecision() const {
	return epsilon;
}

void AugConjugateGradient::setMaximumSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

unsigned int AugConjugateGradient::getMaximumSteps() const {
	return maxSteps;
}

double AugConjugateGradient::getLastError() const {
	return lastError;
}

unsigned int AugConjugateGradient::getLastSteps() const {
	return lastSteps;
}

} /* namespace Update */
