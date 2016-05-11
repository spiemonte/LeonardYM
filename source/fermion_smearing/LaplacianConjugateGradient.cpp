#include "LaplacianConjugateGradient.h"

namespace Update {

LaplacianConjugateGradient::LaplacianConjugateGradient() : epsilon(0.00000000001), maxSteps(30000) { }

LaplacianConjugateGradient::~LaplacianConjugateGradient() { }

bool LaplacianConjugateGradient::solve(const reduced_fermion_lattice_t& lattice, const reduced_color_vector_t& source, reduced_color_vector_t& solution, const block& s, int j_decay, const real_t& shift) {
	solution = source;
	r = zero;
	p = zero;

	Laplacian laplacian(shift, j_decay);
	laplacian.setLattice(lattice);

	reduced_fermion_lattice_t tmp;
	AlgebraUtils::setToZero(tmp);

	laplacian.apply(solution, tmp);

	r = source - tmp;
	p = r;

	real_t norm = AlgebraUtils::squaredNorm(r);
	
	real_t norm_next = norm;

	for (unsigned int step = 0; step < maxSteps; ++step) {
		laplacian.apply(p, tmp);
		
		norm = norm_next;
		std::complex<real_t> alpha = norm/innerProduct(p,tmp);

		solution = solution + alpha*p;
		r = r - alpha*tmp;


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

