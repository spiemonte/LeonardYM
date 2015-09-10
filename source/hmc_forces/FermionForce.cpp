/*
 * FermionForce.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#include "FermionForce.h"
#include "utils/ConvertLattice.h"
#include "utils/ToString.h"

namespace Update {

FermionForce::FermionForce(real_t _kappa) : kappa(_kappa), expMap() { }

FermionForce::~FermionForce() { }

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
		result += -I*imag(trace(derivative*fermionLieGenerator.get(i)*link))*gaugeLieGenerator.get(i);
	}
	return result;
#endif
#ifndef ADJOINT
	GaugeGroup result;
	set_to_zero(result);
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		result += -I*imag(trace(derivative*fermionLieGenerator.get(i)*env.gaugeLinkConfiguration[site][mu]))*gaugeLieGenerator.get(i);
	}
	return result;
#endif
}

/*
ForceVector FermionForce::linkActionDerivative(const environment_t& env, const FermionicForceMatrix& derivative, int site, unsigned int mu) {
	ForceVector result;
	ForceVector omega = expMap.parameters(env.gaugeLinkConfiguration[site][mu]);
	FermionicForceMatrix X;
	set_to_zero(X);
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		X += std::complex<real_t>(0,omega[i])*fermionLieGenerator.get(i);
	}
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		FermionicForceMatrix it = std::complex<real_t>(0,1.)*fermionLieGenerator.get(i), XL = it, XP = X;
		result[i] = 0;
		int factorial = 1;
		for (int n = 0; n < 13; ++n) {
			result[i] += real(trace(derivative*XL))/(factorial);
			XL = it*XP+X*XL;
			XP = XP*X;
			factorial = factorial*(n+2);
		}
	}
	return result;
}

ForceVector FermionForce::linkLieDerivative(const environment_t& env, int site, unsigned int mu, int color) {
	ForceVector result;
	//For every generator
	for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
		result[i] += -2.*real(trace(gaugeLieGenerator.get(i)*gaugeLieGenerator.get(color)*env.gaugeLinkConfiguration[site][mu]));
	}
	return result;
}

GaugeGroup FermionForce::force(const ForceVector& parameters) {
	GaugeGroup result;
	set_to_zero(result);
	//For every generator
	for (unsigned int i = 0; i < gaugeLieGenerator.numberGenerators(); ++i) {
		result += std::complex<real_t>(0,parameters[i])*gaugeLieGenerator.get(i);
	}
	return result;
}*/

} /* namespace Update */
