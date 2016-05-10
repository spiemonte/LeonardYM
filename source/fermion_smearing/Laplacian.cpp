/*
 * Plaquette.cpp
 *
 *  Created on: Feb 29, 2012
 *      Author: spiem_01
 */

#include "Laplacian.h"

namespace Update {

Laplacian::Laplacian(const real_t& _mass, int _j_decay) : mass_sq(_mass), j_decay(_j_decay) { }

Laplacian::~Laplacian() { }

void Laplacian::apply(reduced_color_vector_t& output, const reduced_color_vector_t& input) {
	typedef reduced_fermion_lattice_t LT;
	real_t idf;
	if( j_decay < 4) ftmp = 2.*4-2. + mass_sq;
	else ftmp = 2.*4 + mass_sq;

#pragma omp parallel for
	for (int site = 0; site < input.localsize; ++site) {
		output[site] = ftmp*input[site];
		for (unsigned int mu = 0; mu < 4; ++mu) {
			if (mu != j_decay) output[site] -= environment.gaugeLinkConfiguration[site][mu]*input[LT::sup(site,mu)] + environment.gaugeLinkConfiguration[LT::sdn(site,mu)][mu]*input[LT::sdn(site,mu)];
		}
	}

	output.updateHalo();
}

void Laplacian::setLattice(const reduced_fermion_lattice_t& _lattice) {
	lattice = _lattice;
}

void Laplacian::setMass(const real_t& _mass) {
	mass_sq = _mass;
}

real_t Laplacian::getMass() const {
	return mass_sq;
}

} /* namespace Update */

