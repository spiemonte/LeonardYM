/*
 * Plaquette.h
 *
 *  Created on: Feb 29, 2012
 *      Author: spiem_01
 */

#ifndef LAPLACIAN_H_
#define LAPLACIAN_H_

#include "LatticeSweep.h"

namespace Update {

class Laplacian {
public:
	Laplacian(const real_t& _mass = 0, int _j_decay = 3);
	~Laplacian();

	void apply(reduced_color_vector_t& output, const reduced_color_vector_t& input);

	void setLattice(const reduced_fermion_lattice_t& _lattice);

	void setMass(const real_t& mass);
	real_t getMass() const;
private:
	reduced_fermion_lattice_t lattice;
	real_t mass;
	int j_decay;
};

} /* namespace Update */
#endif /* PLAQUETTE_H_ */

