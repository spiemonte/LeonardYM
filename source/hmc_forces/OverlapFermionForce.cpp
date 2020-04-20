#include "OverlapFermionForce.h"
#include "utils/Gamma.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

OverlapFermionForce::OverlapFermionForce(real_t _kappa, real_t _mass, Polynomial const* _squareRootApproximation) : DiracWilsonFermionForce(_kappa), mass(_mass), squareRootApproximation(_squareRootApproximation) {
	left_dirac_vectors.resize(squareRootApproximation->getRoots().size()+1);
	right_dirac_vectors.resize(squareRootApproximation->getRoots().size()+1);
	diracWilsonOperator = new DiracWilsonOperator();
	squareDiracWilsonOperator = new SquareDiracWilsonOperator();
}

OverlapFermionForce::~OverlapFermionForce() {
	delete diracWilsonOperator;
	delete squareDiracWilsonOperator;
}

void OverlapFermionForce::derivative(extended_fermion_force_lattice_t& fermionForce, const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, real_t weight) {
	diracWilsonOperator->setKappa(kappa);
	squareDiracWilsonOperator->setKappa(kappa);

	diracWilsonOperator->setLattice(lattice);
	squareDiracWilsonOperator->setLattice(lattice);

	std::vector< std::complex<real_t> > roots = squareRootApproximation->getRoots();
	real_t scaling = real(squareRootApproximation->getScaling());

	real_t factor = weight * (0.5-mass/2.);

	//Fermion force is assumed to be initialized outside this routine
	right_dirac_vectors[0] = Y;
	for (unsigned int i = 0; i < roots.size(); ++i) {
		squareDiracWilsonOperator->multiplyAdd(right_dirac_vectors[i+1], right_dirac_vectors[i], right_dirac_vectors[i], -roots[i]);

#pragma omp parallel for
		for (int site = 0; site < right_dirac_vectors[i+1].completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				right_dirac_vectors[i+1][site][mu] = scaling*right_dirac_vectors[i+1][site][mu];
			}
		}
	}

	diracWilsonOperator->multiply(left_dirac_vectors[0], X);
	for (unsigned int i = 0; i < roots.size(); ++i) {
		squareDiracWilsonOperator->multiplyAdd(left_dirac_vectors[i+1], left_dirac_vectors[i], left_dirac_vectors[i], -conj(roots[roots.size() - i - 1]));

#pragma omp parallel for
		for (int site = 0; site < right_dirac_vectors[i+1].completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				left_dirac_vectors[i+1][site][mu] = scaling*left_dirac_vectors[i+1][site][mu];
			}
		}
	}

	for (unsigned int i = 0; i < roots.size(); ++i) {
		diracWilsonOperator->multiply(tmp, left_dirac_vectors[roots.size()-i-1]);

#pragma omp parallel for 
		for (int site = 0; site < fermionForce.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				fermionForce[site][mu] -= scaling * factor * (DiracWilsonFermionForce::derivative(lattice, tmp, right_dirac_vectors[i], site, mu));
			}
		}

		diracWilsonOperator->multiply(tmp, right_dirac_vectors[i]);

#pragma omp parallel for 
		for (int site = 0; site < fermionForce.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				fermionForce[site][mu] -= scaling * factor * (DiracWilsonFermionForce::derivative(lattice, left_dirac_vectors[roots.size()-i-1], tmp, site, mu));
			}
		}
	}

#pragma omp parallel for 
	for (int site = 0; site < fermionForce.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			fermionForce[site][mu] -= factor * (DiracWilsonFermionForce::derivative(lattice, X, right_dirac_vectors.back(), site, mu));
		}
	}
}

} /* namespace Update */
