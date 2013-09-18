/*
 * BlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "BlockDiracWilsonOperator.h"
#include "../BlockDiracWilsonFermionForce.h"

namespace Update {

BlockDiracWilsonOperator::BlockDiracWilsonOperator() : BlockDiracOperator()  { }

BlockDiracWilsonOperator::BlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa) : BlockDiracOperator(_lattice, _kappa) { }

BlockDiracWilsonOperator::~BlockDiracWilsonOperator() { }

void BlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	typedef reduced_fermion_lattice_t LT;
	typedef reduced_dirac_vector_t DV;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (unsigned int site = 0; site < lattice.localsize; ++site) {
		GaugeVector tmp[4][4];

		//First we put U(x,mu)*input(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = lattice[site][mu]*input[DV::sup(site,mu)][nu];
			}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4];
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] = (tmp[0][0][n]+I*tmp[0][3][n]+tmp[1][0][n]+tmp[1][3][n]+tmp[2][0][n]+I*tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[1][n] = (tmp[0][1][n]+I*tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+tmp[2][1][n]-I*tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
			hopping[2][n] = (I*tmp[0][1][n]-tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+I*tmp[2][0][n]-tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] = (I*tmp[0][0][n]-tmp[0][3][n]-tmp[1][0][n]-tmp[1][3][n]-I*tmp[2][1][n]-tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
		}

		//Then we put U(x-mu,mu)*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = htrans(lattice[LT::sdn(site,mu)][mu])*input[DV::sdn(site,mu)][nu];
			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] += (tmp[0][0][n]-I*tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+tmp[2][0][n]-I*tmp[2][2][n]+tmp[3][0][n]+tmp[3][2][n]);
			hopping[1][n] += (tmp[0][1][n]-I*tmp[0][2][n]+tmp[1][1][n]+tmp[1][2][n]+tmp[2][1][n]+I*tmp[2][3][n]+tmp[3][1][n]+tmp[3][3][n]);
			hopping[2][n] += (-I*tmp[0][1][n]-tmp[0][2][n]-tmp[1][1][n]-tmp[1][2][n]-I*tmp[2][0][n]-tmp[2][2][n]-tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] += (-I*tmp[0][0][n]-tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+I*tmp[2][1][n]-tmp[2][3][n]-tmp[3][1][n]-tmp[3][3][n]);
		}

		//The final result is gamma5*input - kappa*hopping
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			output[site][0][n] = + input[site][0][n] - kappa*hopping[0][n];
			output[site][1][n] = + input[site][1][n] - kappa*hopping[1][n];
			output[site][2][n] = - (input[site][2][n] + kappa*hopping[2][n]);
			output[site][3][n] = - (input[site][3][n] + kappa*hopping[3][n]);
		}
	}
	output.updateHalo();
}

void BlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	typedef reduced_fermion_lattice_t LT;
	typedef reduced_dirac_vector_t DV;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (unsigned int site = 0; site < lattice.localsize; ++site) {
		//First we start the hopping parameter terms
		GaugeVector tmp[4][4];

		//First we put U(x,mu)*input(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = lattice[site][mu]*vector1[DV::sup(site,mu)][nu];
			}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4];
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] = (tmp[0][0][n]+I*tmp[0][3][n]+tmp[1][0][n]+tmp[1][3][n]+tmp[2][0][n]+I*tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[1][n] = (tmp[0][1][n]+I*tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+tmp[2][1][n]-I*tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
			hopping[2][n] = (I*tmp[0][1][n]-tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+I*tmp[2][0][n]-tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] = (I*tmp[0][0][n]-tmp[0][3][n]-tmp[1][0][n]-tmp[1][3][n]-I*tmp[2][1][n]-tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
		}

		//Then we put U(x-mu,mu)*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = htrans(lattice[LT::sdn(site,mu)][mu])*vector1[DV::sdn(site,mu)][nu];
			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] += (tmp[0][0][n]-I*tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+tmp[2][0][n]-I*tmp[2][2][n]+tmp[3][0][n]+tmp[3][2][n]);
			hopping[1][n] += (tmp[0][1][n]-I*tmp[0][2][n]+tmp[1][1][n]+tmp[1][2][n]+tmp[2][1][n]+I*tmp[2][3][n]+tmp[3][1][n]+tmp[3][3][n]);
			hopping[2][n] += (-I*tmp[0][1][n]-tmp[0][2][n]-tmp[1][1][n]-tmp[1][2][n]-I*tmp[2][0][n]-tmp[2][2][n]-tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] += (-I*tmp[0][0][n]-tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+I*tmp[2][1][n]-tmp[2][3][n]-tmp[3][1][n]-tmp[3][3][n]);
		}

		//The final result is gamma5*input - kappa*hopping
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			output[site][0][n] = alpha*vector2[site][0][n] + vector1[site][0][n] - kappa*hopping[0][n];
			output[site][1][n] = alpha*vector2[site][1][n] + vector1[site][1][n] - kappa*hopping[1][n];
			output[site][2][n] = alpha*vector2[site][2][n] - (vector1[site][2][n] + kappa*hopping[2][n]);
			output[site][3][n] = alpha*vector2[site][3][n] - (vector1[site][3][n] + kappa*hopping[3][n]);
		}
	}
	output.updateHalo();//TODO is needed?
}

void BlockDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	typedef reduced_fermion_lattice_t::Layout Layout;
	for (int site = 0; site < lattice.localsize; ++site) {
		if (Layout::globalIndexX(site) % blockSize == 0) set_to_zero(this->lattice[site][0]);
		if (Layout::globalIndexY(site) % blockSize == 0) set_to_zero(this->lattice[site][1]);
		if (Layout::globalIndexZ(site) % blockSize == 0) set_to_zero(this->lattice[site][2]);
		if (Layout::globalIndexT(site) % blockSize == 0) set_to_zero(this->lattice[site][3]);
	}
}

FermionForce* BlockDiracWilsonOperator::getForce() const {
	return new BlockDiracWilsonFermionForce(kappa, blockSize);
}

} /* namespace Update */
