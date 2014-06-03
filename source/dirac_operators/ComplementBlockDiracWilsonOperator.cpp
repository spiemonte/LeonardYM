/*
 * ComplementBlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "ComplementBlockDiracWilsonOperator.h"
#include "../DiracWilsonFermionForce.h"
#include "../AlgebraUtils.h"

namespace Update {

ComplementBlockDiracWilsonOperator::ComplementBlockDiracWilsonOperator() : BlockDiracOperator(), diracWilsonOperator(), dir(0) { }

ComplementBlockDiracWilsonOperator::ComplementBlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : BlockDiracOperator(_lattice, _kappa), diracWilsonOperator(_lattice, _kappa), dir(0) {
	this->setLattice(_lattice);
}

ComplementBlockDiracWilsonOperator::~ComplementBlockDiracWilsonOperator() { }

void ComplementBlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	diracWilsonOperator.multiply(output, input);
}

void ComplementBlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	diracWilsonOperator.multiplyAdd(output, vector1, vector2, alpha);
}

void ComplementBlockDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	typedef reduced_fermion_lattice_t::Layout Layout;
	for (int site = 0; site < lattice.localsize; ++site) {
		if (Layout::globalIndexX(site) % xBlockSize != 0 || (Layout::globalIndexX(site) % xBlockSize == 0 && dir != 0) ) set_to_zero(this->lattice[site][0]);
		if (Layout::globalIndexY(site) % yBlockSize != 0 || (Layout::globalIndexY(site) % yBlockSize == 0 && dir != 1) ) set_to_zero(this->lattice[site][1]);
		if (Layout::globalIndexZ(site) % zBlockSize != 0 || (Layout::globalIndexZ(site) % zBlockSize == 0 && dir != 2) ) set_to_zero(this->lattice[site][2]);
		if (Layout::globalIndexT(site) % tBlockSize != 0 || (Layout::globalIndexT(site) % tBlockSize == 0 && dir != 3) ) set_to_zero(this->lattice[site][3]);
	}
	diracWilsonOperator.setLattice(this->lattice);
}

FermionForce* ComplementBlockDiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

void ComplementBlockDiracWilsonOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracWilsonOperator.setKappa(_kappa);
}

void ComplementBlockDiracWilsonOperator::setDirection(unsigned int mu) {
	dir = mu;
}

unsigned int ComplementBlockDiracWilsonOperator::getDirection() const {
	return dir;
}


} /* namespace Update */
