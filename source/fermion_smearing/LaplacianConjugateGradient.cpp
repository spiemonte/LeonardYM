#include "LaplacianConjugateGradient.h"

namespace Update {

LaplacianConjugateGradient::LaplacianConjugateGradient() : epsilon(0.00000000001), maxSteps(30000) { }

LaplacianConjugateGradient::~LaplacianConjugateGradient() { }

bool LaplacianConjugateGradient::solve(const multi1d<LatticeColorMatrix>& u, const LatticeColorVector &source, LatticeColorVector &solution, const Subset & s, int j_decay, const Double& shift) {
	solution = source;
	r = zero;
	p = zero;

	LatticeColorVector tmp = zero;
	klein_gord(u, solution, tmp, shift, j_decay);

	r[s] = source - tmp;
	p[s] = r;

	Double norm = norm2(r);
	
	Double norm_next = norm;

	for (unsigned int step = 0; step < maxSteps; ++step) {
		//tmp[s] = p;
		//laplacian(u, tmp, j_decay, 1);
		klein_gord(u, p, tmp, shift, j_decay);
		
		norm = norm_next;
		Complex alpha = norm/innerProduct(p,tmp);

		solution[s] = solution + alpha*p;
		r[s] = r - alpha*tmp;


		norm_next = norm2(r);
		if (toBool(norm_next < epsilon)) {
			QDPIO::cout << "LaplacianConjugateGradient::Convergence in " << step << " steps, final norm: " << norm_next << std::endl;
			lastSteps = step;
			return true;
		}
		/*else {
			QDPIO::cout << "LaplacianConjugateGradient::Residual norm at step " << step << ": " << norm_next << std::endl;
		}*/

		Double beta = norm_next/norm;
		
		p[s] = r + beta*p;
	}

	QDPIO::cout << "LaplacianConjugateGradient::Failure in finding convergence, last error: " << norm_next << std::endl;
	return false;
}

void LaplacianConjugateGradient::setPrecision(Real _epsilon) {
	epsilon = _epsilon;
}

Real LaplacianConjugateGradient::getPrecision() const {
	return epsilon;
}

void LaplacianConjugateGradient::setMaximumSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

unsigned int LaplacianConjugateGradient::getMaximumSteps() const {
	return maxSteps;
}

Real LaplacianConjugateGradient::getLastError() const {
	return lastError;
}

unsigned int LaplacianConjugateGradient::getLastSteps() const {
	return lastSteps;
}

}

