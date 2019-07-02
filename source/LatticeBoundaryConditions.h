#ifndef LATTICEBOUNDARYCONDITIONS_H
#define LATTICEBOUNDARYCONDITIONS_H

#include "LatticeDefinitions.h"

//Functions to enforce the given boundary conditions for fermions and gauge fields

template<typename T> void switchAntiperiodicBc(T& lattice) {
	typedef typename T::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		if (LT::globalIndexT(site) == (LT::glob_t - 1)) lattice[site][3] = -lattice[site][3];
	}
}

template<typename T> void switchSpatialAntiperiodicBc(T& lattice) {
	typedef typename T::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		if (LT::globalIndexZ(site) == (LT::glob_z - 1)) lattice[site][2] = -lattice[site][2];
	}
}

template<typename T> void switchFullAntiperiodicBc(T& lattice) {
        typedef typename T::Layout LT;
#pragma omp parallel for
        for (int site = 0; site < lattice.completesize; ++site) {
		if (LT::globalIndexX(site) == (LT::glob_x - 1)) lattice[site][0] = -lattice[site][0];
		if (LT::globalIndexY(site) == (LT::glob_y - 1)) lattice[site][1] = -lattice[site][1];
		if (LT::globalIndexZ(site) == (LT::glob_z - 1)) lattice[site][2] = -lattice[site][2];
		if (LT::globalIndexT(site) == (LT::glob_t - 1)) lattice[site][3] = -lattice[site][3];
        }
}

template<typename T> void switchOpenBc(T& lattice) {
	typedef typename T::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		int t_index = LT::globalIndexT(site);
		if (t_index == (LT::glob_t - 1)) {
			for (unsigned int mu = 0; mu < 4; ++mu) set_to_zero(lattice[site][mu]);
		}
		else if (t_index == (LT::glob_t - 2)) {
			set_to_zero(lattice[site][3]);
		}
	}
}

#ifdef ADJOINT
inline void switchOpenBc(extended_fermion_lattice_t& lattice) {
	typedef extended_fermion_lattice_t::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		int t_index = LT::globalIndexT(site);
		if (t_index == (LT::glob_t - 1)) {
			for (unsigned int mu = 0; mu < 4; ++mu) set_to_zero(lattice[site][mu]);
		}
		else if (t_index == (LT::glob_t - 2)) {
			set_to_zero(lattice[site][3]);
		}
	}
}
#endif

#endif
