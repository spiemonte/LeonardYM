#ifndef DIRACWILSONFERMIONFORCE_H_
#define DIRACWILSONFERMIONFORCE_H_

#include "FermionForce.h"

namespace Update {

class DiracWilsonFermionForce: public Update::FermionForce {
public:
	DiracWilsonFermionForce(real_t _kappa);
	~DiracWilsonFermionForce();

	virtual FermionicForceMatrix derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const;
};

} /* namespace Update */
#endif /* DIRACWILSONFERMIONFORCE_H_ */
