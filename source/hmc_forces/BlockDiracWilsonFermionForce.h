/*
 * BlockDiracWilsonFermionForce.h
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#ifndef BLOCKDIRACWILSONFERMIONFORCE_H_
#define BLOCKDIRACWILSONFERMIONFORCE_H_

#include "FermionForce.h"

namespace Update {

class BlockDiracWilsonFermionForce: public Update::FermionForce {
	public:
		BlockDiracWilsonFermionForce(real_t _kappa, int _xBlockSize, int _yBlockSize, int _zBlockSize, int _tBlockSize);
		~BlockDiracWilsonFermionForce();

		virtual FermionicForceMatrix derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const;
	private:
		int xBlockSize;
		int yBlockSize;
		int zBlockSize;
		int tBlockSize;
};

} /* namespace Update */

#endif
