// Author: Stefano Piemonte, (C) 2017
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef FERMIONFORCE_H_
#define FERMIONFORCE_H_
#include "Environment.h"
#include "utils/LieGenerators.h"
#include "utils/ExpMap.h"

namespace Update {

template <class TLinkConfig, class TFermionLinkConfig, class TVector> class FermionForce {
public:
	typedef typename TLinkConfig::TMatrix TMatrix;
	typedef typename TLinkConfig::MatrixLib MatrixLib;
	typedef typename TLinkConfig::Real Real;
	//We need full color complex matrices even in the adjoint case
	typedef typename PlainMatrixLibrary<std::complex<Real>, TVector::grouplength>::TMatrix FermionicForceMatrix;
	typedef typename PlainMatrixLibrary<std::complex<Real>, TVector::grouplength> FermionicForceMatrixLib;
	//Color vector
	typedef typename TVector::GroupVector TGroupVector;

	FermionForce(Real _kappa) : kappa(_kappa) { }
	virtual ~FermionForce() { }

	virtual FermionicForceMatrix derivative(const TFermionLinkConfig& lattice, const TVector& X, const TVector& Y, int site, int mu) const = 0;

	TMatrix force(const TFermionLinkConfig& fermionLinkConfig, const FermionicForceMatrix& derivative, int site, unsigned int mu) {
		TMatrix result;
		MatrixLib::setToZero(result);
		//For every generator
		for (unsigned int i = 0; i < fermionLieGenerator.numberGenerators(); ++i) {
			result += -I*imag(trace(derivative*fermionLieGenerator.get(i)*fermionLinkConfig[site][mu]))*gaugeLieGenerator.get(i);
		}
		return result;
	}

	//Required by the clover term force
	virtual void setLattice(const TFermionLinkConfig& ) { }

protected:
	Real kappa;
	
	//ExponentialMap expMap;

	FermionicForceMatrix tensor(const TGroupVector& x, const TGroupVector& y) const {
		FermionicForceMatrix result;
		FermionicForceMatrixLib::setToZero(result);
		for (int i = 0; i < TVector::grouplength; ++i) {
			for (int j = 0; j < TVector::grouplength; ++j) {
				result(i,j) = y[i]*conj(x[j]);
			}
		}
		return result;
	}

	LieGenerator<TFermionLinkConfig::TMatrix> fermionLieGenerator;
	LieGenerator<TLinkConfig::TMatrix> gaugeLieGenerator;
};

} /* namespace Update */
#endif /* FERMIONFORCE_H_ */
