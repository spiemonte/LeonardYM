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
#include "./MPILattice/ReducedUpStencil.h"
#include "./MPILattice/ReducedDownStencil.h"
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

typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::MpiLayout<Lattice::ReducedUpStencil> > half_spinor_up_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::MpiLayout<Lattice::ReducedDownStencil> > half_spinor_down_vector_t;

typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::StandardStencil> > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_field_strength_lattice_t;

typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::StandardStencil> > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_matrix_lattice_t;

typedef Lattice::Lattice<int[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_index_lattice_t;

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

typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_extended_dirac_vector_t;
typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_standard_dirac_vector_t;
typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_reduced_dirac_vector_t;

typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::LocalLayout > half_spinor_up_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::LocalLayout > half_spinor_down_vector_t;

typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > reduced_field_strength_lattice_t;

typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > reduced_matrix_lattice_t;

typedef Lattice::Lattice<int[4], Lattice::LocalLayout > extended_index_lattice_t;

#endif

#ifdef ENABLE_MPI
typedef Lattice::Lattice<int, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_index_lattice_t;
#endif
#ifndef ENABLE_MPI
typedef Lattice::Lattice<int, Lattice::LocalLayout > reduced_index_lattice_t;
#endif

inline void switchAntiperiodicBc(extended_fermion_lattice_t& lattice) {
	typedef extended_fermion_lattice_t::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		if (LT::globalIndexT(site) == (LT::glob_t - 1)) lattice[site][3] = -lattice[site][3];
	}
	//lattice.updateHalo(); TODO
}

inline void switchSpatialAntiperiodicBc(extended_fermion_lattice_t& lattice) {
	typedef extended_fermion_lattice_t::Layout LT;
#pragma omp parallel for
	for (int site = 0; site < lattice.completesize; ++site) {
		if (LT::globalIndexZ(site) == (LT::glob_z - 1)) lattice[site][2] = -lattice[site][2];
	}
	//lattice.updateHalo(); TODO
}

inline void switchOpenBc(extended_gauge_lattice_t& lattice) {
	typedef extended_gauge_lattice_t::Layout LT;
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
	//lattice.updateHalo();
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
	//lattice.updateHalo();
}
#endif

#include "StorageParameters.h"
#include "ConvertLattice.h"

namespace Update {

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
		try {
			std::string bc = configurations.get<std::string>("boundary_conditions");
			if (bc == "fermion_antiperiodic") {
				switchAntiperiodicBc(fermionicLinkConfiguration);
			}
			else if (bc == "fermion_spatial_antiperiodic") {
				switchSpatialAntiperiodicBc(fermionicLinkConfiguration);
			}
			else if (bc == "fermion_periodic") {

			}
			else if (bc == "open") {
				switchOpenBc(fermionicLinkConfiguration);
				switchOpenBc(gaugeLinkConfiguration);
			}
			else {
				static int count = 0;
				if (count == 0 && isOutputProcess()) {
					std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
					count = count + 1;
				}
				switchAntiperiodicBc(fermionicLinkConfiguration);
			}
		}
		catch (NotFoundOption& ex) {
			static int count = 0;
			if (count == 0 && isOutputProcess()) {
				std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
				count = count + 1;
			}
			switchAntiperiodicBc(fermionicLinkConfiguration);
		}
	}

	void setFermionBc(extended_fermion_lattice_t& lattice) const {
		try {
			std::string bc = configurations.get<std::string>("boundary_conditions");
			if (bc == "fermion_antiperiodic") {
				switchAntiperiodicBc(lattice);
			}
			else if (bc == "fermion_spatial_antiperiodic") {
				switchSpatialAntiperiodicBc(lattice);
			}
			else if (bc == "fermion_periodic") {

			}
			else if (bc == "open") {
				switchOpenBc(lattice);
			}
			else {
				static int count = 0;
				if (count == 0 && isOutputProcess()) {
					std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
					count = count + 1;
				}
				switchAntiperiodicBc(lattice);
			}
		}
		catch (NotFoundOption& ex) {
			static int count = 0;
			if (count == 0 && isOutputProcess()) {
				std::cout << "Warning, boundary conditions not or badly set, using antiperiodic!" << std::endl;
				count = count + 1;
			}
			switchAntiperiodicBc(lattice);
		}
	}

	Environment& operator=(const Environment& toCopy) {
		configurations = toCopy.configurations;
		sweep = toCopy.sweep;
		iteration = toCopy.iteration;
		measurement = toCopy.measurement;
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
inline void reduceAllSum(float& value) {
	float result = value;
	MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(float&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(int& value) {
	int result = value;
	MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(int&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(std::complex<double>& value) {
	std::complex<double> result = value;
	MPI_Allreduce(&value.real(), &result.real(), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&value.imag(), &result.imag(), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(std::complex<double>&) { }
#endif

#ifdef ENABLE_MPI
inline void reduceAllSum(std::complex<float>& value) {
	std::complex<float> result = value;
	MPI_Allreduce(&value.real(), &result.real(), 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&value.imag(), &result.imag(), 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
	value = result;
}
#endif
#ifndef ENABLE_MPI
inline void reduceAllSum(std::complex<float>&) { }
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
