/*
 * DiracWilsonFermionForce.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: spiem_01
 */

#include "BlockDiracWilsonFermionForce.h"
#include "Gamma.h"
#include "ToString.h"

namespace Update {

BlockDiracWilsonFermionForce::BlockDiracWilsonFermionForce(real_t _kappa, int _xBlockSize, int _yBlockSize, int _zBlockSize, int _tBlockSize) : FermionForce(_kappa), xBlockSize(_xBlockSize), yBlockSize(_yBlockSize), zBlockSize(_zBlockSize), tBlockSize(_tBlockSize) { }

BlockDiracWilsonFermionForce::~BlockDiracWilsonFermionForce() { }

FermionicForceMatrix BlockDiracWilsonFermionForce::derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const {
	typedef extended_dirac_vector_t Vector;
	typedef extended_fermion_lattice_t::Layout Layout;
	FermionicForceMatrix force;
	set_to_zero(force);

	//We have to ensure that we are not in the block boundary
	bool inBlockBoundary = false;
	if (mu == 0 && Layout::globalIndexX(site) % xBlockSize == 0) inBlockBoundary = true;
	if (mu == 1 && Layout::globalIndexY(site) % yBlockSize == 0) inBlockBoundary = true;
	if (mu == 2 && Layout::globalIndexZ(site) % zBlockSize == 0) inBlockBoundary = true;
	if (mu == 3 && Layout::globalIndexT(site) % tBlockSize == 0) inBlockBoundary = true;

	if (inBlockBoundary) {
		return force;
	}
	else {
		//Derivative from X(x)^\dag\gamma_5(1-\gamma_\mu) U_\mu(x)Y(x+mu)
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			GaugeVector tmp;
			set_to_zero(tmp);
			for (unsigned int beta = 0; beta < 4; ++beta) {
				if (Gamma::g5idmg(mu,alpha,beta) != static_cast<real_t>(0.)) tmp += -(kappa)*Gamma::g5idmg(mu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
			}
			force += tensor(X[site][alpha],tmp);
		}
		//Derivative from X(x+mu)^\dag\gamma_5(1+\gamma_\mu) U^\dag_\mu(x)Y(x)
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			GaugeVector tmp;
			set_to_zero(tmp);
			for (unsigned int beta = 0; beta < 4; ++beta) {
				if (Gamma::g5idpg(mu,alpha,beta) != static_cast<real_t>(0.)) tmp += -(kappa)*Gamma::g5idpg(mu,alpha,beta)*htrans(lattice[site][mu])*Y[site][beta];
			}
			force -= tensor(lattice[site][mu]*X[Vector::sup(site,mu)][alpha],tmp);
		}

		//Now we change just X with Y
		//Derivative from Y(x)^\dag\gamma_5(1-\gamma_\mu) U_\mu(x)X(x+mu)
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			GaugeVector tmp;
			set_to_zero(tmp);
			for (unsigned int beta = 0; beta < 4; ++beta) {
				if (Gamma::g5idmg(mu,alpha,beta) != static_cast<real_t>(0.)) tmp += -(kappa)*Gamma::g5idmg(mu,alpha,beta)*X[Vector::sup(site,mu)][beta];
			}
			force += tensor(Y[site][alpha],tmp);
		}
		//Derivative from Y(x+mu)^\dag\gamma_5(1+\gamma_\mu) U^\dag_\mu(x)X(x)
		for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			GaugeVector tmp;
			set_to_zero(tmp);
			for (unsigned int beta = 0; beta < 4; ++beta) {
				if (Gamma::g5idpg(mu,alpha,beta) != static_cast<real_t>(0.)) tmp += -(kappa)*Gamma::g5idpg(mu,alpha,beta)*htrans(lattice[site][mu])*X[site][beta];
			}
			force -= tensor(lattice[site][mu]*Y[Vector::sup(site,mu)][alpha],tmp);
		}

		return force;
	}
}

} /* namespace Update */
