/*
 * Enviroment.h
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include "MatrixTypedef.h"
#include "./MPILattice/Lattice.h"

#ifdef ENABLE_MPI
#include "./MPILattice/MPILayout.h"
#include "./MPILattice/ReducedStencil.h"
#include "./MPILattice/StandardStencil.h"
#include "./MPILattice/ExtendedStencil.h"

typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_gauge_lattice_t;

typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_fermion_lattice_t;

typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_dirac_vector_t;

typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::StandardStencil> > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_field_strength_lattice_t;

typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::StandardStencil> > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_matrix_lattice_t;

#endif
#ifndef ENABLE_MPI
#include "./MPILattice/LocalLayout.h"

typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > extended_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > standard_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > reduced_gauge_lattice_t;

typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > extended_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > standard_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > reduced_fermion_lattice_t;

typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > extended_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > standard_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > reduced_dirac_vector_t;

typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > reduced_field_strength_lattice_t;

typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > reduced_matrix_lattice_t;

#endif

inline void switchAntiperiodicBc(extended_fermion_lattice_t& lattice) {
	typedef extended_fermion_lattice_t::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		if (LT::globalIndexT(site) == (LT::glob_t - 1)) lattice[site][3] = -lattice[site][3];
	}
	lattice.updateHalo();
}

#include "StorageParameters.h"
#include "ConvertLattice.h"

namespace Update {

// @class Environment
// @brief This class is a container of the dynamical

class Environment {
public:
	Environment() : gaugeLinkConfiguration(), fermionicLinkConfiguration(), sweep(0), iteration(0), measurement(false) { }
	Environment(const boost::program_options::variables_map& vm) : gaugeLinkConfiguration(), fermionicLinkConfiguration(), configurations(vm), sweep(0), iteration(0), measurement(false) { }
	Environment(const Environment& toCopy) : gaugeLinkConfiguration(toCopy.gaugeLinkConfiguration), fermionicLinkConfiguration(toCopy.fermionicLinkConfiguration), configurations(toCopy.configurations), sweep(toCopy.sweep), iteration(toCopy.iteration), measurement(toCopy.measurement) { }
	~Environment() { }

	void synchronize() {
#ifdef ADJOINT
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::convert(fermionicLinkConfiguration, gaugeLinkConfiguration);
#endif
#ifndef ADJOINT
		fermionicLinkConfiguration = gaugeLinkConfiguration;
#endif
		switchAntiperiodicBc(fermionicLinkConfiguration);
	}

	Environment& operator=(const Environment& toCopy) {
		gaugeLinkConfiguration = toCopy.gaugeLinkConfiguration;
		fermionicLinkConfiguration = toCopy.fermionicLinkConfiguration;
		return *this;
	}

	extended_fermion_lattice_t& getFermionLattice() {
		return fermionicLinkConfiguration;
	}

	const extended_fermion_lattice_t& getFermionLattice() const {
		return fermionicLinkConfiguration;
	}

	extended_gauge_lattice_t gaugeLinkConfiguration;
	extended_fermion_lattice_t fermionicLinkConfiguration;
	StorageParameters configurations;

	unsigned int sweep;
	unsigned int iteration;

	bool measurement;
};

inline bool isOutputProcess() {
#ifdef ENABLE_MPI
	int this_processor;
	MPI_Comm_rank(MPI_COMM_WORLD, &this_processor);
	return this_processor == 0;
#endif
#ifndef ENABLE_MPI
	return true;
#endif
}

#ifdef ENABLE_MPI
inline void reduceAllSum(double& value) {
	double result = value;
	MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(double&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(long double& value) {
	long double result = value;
	MPI_Allreduce(&value, &result, 1, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(long double&) { }
#endif


typedef Environment environment_t;

}

#endif /* ENVIROMENT_H_ */
