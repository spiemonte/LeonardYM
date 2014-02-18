/*
 * ImprovedFermionForce.h
 *
 *  Created on: Jun 11, 2012
 *      Author: spiem_01
 */

#ifndef IMPROVEDFERMIONFORCE_H_
#define IMPROVEDFERMIONFORCE_H_

#include "DiracWilsonFermionForce.h"

namespace Update {

//This class stores the part depending on the lattice for a single term of clover force
class CloverLatticeForce {
public:
	CloverLatticeForce() : memory(0), size(0) {	}

	CloverLatticeForce(const CloverLatticeForce& snd) {
		if (snd.memory != 0) {
			size = snd.size;
			memory = new LinkForceElement[size];
#pragma omp parallel for
			for (int site = 0; site < size; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					for (unsigned int nu = 0; nu < 4; ++nu) {
						memory[site][mu][nu] = snd.memory[site][mu][nu];
					}
				}
			}
		}
		else {
			memory = 0;
			size = 0;
		}
	}

	void setLocaLatticeSize(int _size) {
		if (size != _size) {
			size = _size;
			if (memory != 0) delete[] memory;
			memory = new LinkForceElement[size];
		}
	}

	~CloverLatticeForce() {
		if (memory != 0) delete[] memory;
	}

	FermionicGroup& operator()(int site, int mu, int nu) {
		return memory[site][mu][nu];
	}

	const FermionicGroup& operator()(int site, int mu, int nu) const {
		return memory[site][mu][nu];
	}

private:
	typedef FermionicGroup LinkForceElement[4][4];

	LinkForceElement* memory;
	int size;
};

class ImprovedFermionForce : public DiracWilsonFermionForce {
public:
	ImprovedFermionForce(real_t _kappa, real_t _csw = 1.);
	ImprovedFermionForce(const ImprovedFermionForce& toCopy);
	~ImprovedFermionForce();

	virtual FermionicForceMatrix derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const;

	virtual void setLattice(const extended_fermion_lattice_t& lattice);
private:
	real_t csw;

	CloverLatticeForce latticeForcesLeft[8];
	CloverLatticeForce latticeForcesRight[8];
};

} /* namespace Update */
#endif /* IMPROVEDFERMIONFORCE_H_ */
