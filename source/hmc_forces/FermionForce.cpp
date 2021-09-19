#include "FermionForce.h"
#include "utils/ConvertLattice.h"
#include "utils/ToString.h"

namespace Update {

FermionForce::FermionForce(real_t _kappa) : kappa(_kappa), expMap() { }

FermionForce::~FermionForce() { }

void FermionForce::derivative(extended_fermion_force_lattice_t& fermionForce, const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, real_t weight) {
#pragma omp parallel for
	for (int site = 0; site < fermionForce.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			//Minus sign on the fermion force!
			fermionForce[site][mu] -= weight * (this->derivative(lattice, X, Y, site, mu));
		}
	}

}

FermionicForceMatrix FermionForce::tensor(const GaugeVector& x, const GaugeVector& y) const {
	FermionicForceMatrix result;
	set_to_zero(result);
	for (int i = 0; i < diracVectorLength; ++i) {
		for (int j = 0; j < diracVectorLength; ++j) {
			result.at(i,j) = y[i]*conj(x[j]);
		}
	}
	return result;
}

GaugeGroup FermionForce::force(const environment_t& env, const FermionicForceMatrix& derivative, int site, unsigned int mu) {
#ifdef ADJOINT
	GaugeGroup result;
	set_to_zero(result);
	FermionicGroup link;
	ConvertLattice<GaugeGroup,FermionicGroup>::toAdjoint(env.gaugeLinkConfiguration[site][mu],link);
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		result += -II*imag(trace(derivative*fermionLieGenerator.get(i)*link))*gaugeLieGenerator.get(i);
	}
	return result;
#endif
#ifndef ADJOINT
	GaugeGroup result;
	set_to_zero(result);
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		result += -II*imag(trace(derivative*fermionLieGenerator.get(i)*env.gaugeLinkConfiguration[site][mu]))*gaugeLieGenerator.get(i);
	}
	return result;
#endif
}

} /* namespace Update */
