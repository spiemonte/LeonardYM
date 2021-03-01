#ifndef LATTICEDEFINITIONS_H
#define LATTICEDEFINITIONS_H

#include "MatrixTypedef.h"
#include "./MPILattice/Lattice.h"

#ifdef ALIGNED_OPT
#include "./MPILattice/MDLattice.h"
#endif

//Typedefs of all lattice fields required by the Monte Carlo simulations of QCD-like theories

#ifdef ENABLE_MPI
#include "./MPILattice/MPILayout.h"
#include "./MPILattice/ReducedStencil.h"
#include "./MPILattice/ReducedUpStencil.h"
#include "./MPILattice/ReducedDownStencil.h"
#include "./MPILattice/StandardStencil.h"
#include "./MPILattice/ExtendedStencil.h"

//Gauge field configurations for different MPI layouts
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_gauge_lattice_t;

//Force field configurations for different MPI layouts
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_force_lattice_t;
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_force_lattice_t;
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_force_lattice_t;

//Gauge field configurations interacting with fermions for different MPI layouts
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_fermion_lattice_t;

//Fermion force field
typedef Lattice::Lattice<Update::FermionicForceMatrix[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_fermion_force_lattice_t;

//(Pseudo)Fermion field configurations for different MPI layouts
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_dirac_vector_t;

//Half spinors projections for reduced layouts
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::MpiLayout<Lattice::ReducedUpStencil> > half_spinor_up_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::MpiLayout<Lattice::ReducedDownStencil> > half_spinor_down_vector_t;

//Gauge field strength F_\mu\nu for different MPI layouts
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::StandardStencil> > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_field_strength_lattice_t;

//Gauge transformation lattice for different MPI layouts
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::StandardStencil> > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_matrix_lattice_t;

//Gauge field configuration in the adjoint representations for different MPI layouts
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_adjoint_lattice_t;
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::MpiLayout<Lattice::StandardStencil> > standard_adjoint_lattice_t;
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_adjoint_lattice_t;

//Scalar field in the fundamental representation of the gauge group for different MPI layouts
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_color_vector_t;
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::MpiLayout<Lattice::StandardStencil> > standard_color_vector_t;
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_color_vector_t;

//Scalar field in the adjoint representation of the gauge group for different MPI layouts
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_adjoint_complex_color_vector_t;
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::MpiLayout<Lattice::StandardStencil> > standard_adjoint_complex_color_vector_t;
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_adjoint_complex_color_vector_t;

//Scalar field in the adjoint representation of the gauge group for different MPI layouts
typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_adjoint_real_color_vector_t;
typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::MpiLayout<Lattice::StandardStencil> > standard_adjoint_real_color_vector_t;
typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_adjoint_real_color_vector_t;

//Lattice of indeces to store positions and similar variables
typedef Lattice::Lattice<int[4], Lattice::MpiLayout<Lattice::ExtendedStencil> > extended_index_lattice_t;

#ifdef ALIGNED_OPT
typedef Lattice::ComplexVectorLattice<Update::real_t, Update::GaugeVector[4], Lattice::MpiLayout<Lattice::ReducedStencil>, 4, Update::diracVectorLength> reduced_soa_dirac_vector_t;
#endif

#endif
#ifndef ENABLE_MPI
#include "./MPILattice/LocalLayout.h"

//Gauge field configurations equal for all layouts
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > extended_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > standard_gauge_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup[4], Lattice::LocalLayout > reduced_gauge_lattice_t;

//Force field configurations equal for all layouts
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::LocalLayout > extended_force_lattice_t;
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::LocalLayout > standard_force_lattice_t;
typedef Lattice::Lattice<Update::ForceVector[4], Lattice::LocalLayout > reduced_force_lattice_t;

//Gauge field configurations interacting with fermions equal for all layouts
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > extended_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > standard_fermion_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[4], Lattice::LocalLayout > reduced_fermion_lattice_t;

//Fermion force field
typedef Lattice::Lattice<Update::FermionicForceMatrix[4], Lattice::LocalLayout > extended_fermion_force_lattice_t;

//(Pseudo)Fermion field configurations equal for all layouts
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > extended_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > standard_dirac_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[4], Lattice::LocalLayout > reduced_dirac_vector_t;

typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_extended_dirac_vector_t;
typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_standard_dirac_vector_t;
typedef Lattice::Lattice<Update::single_GaugeVector[4], Lattice::LocalLayout > single_reduced_dirac_vector_t;

//Half spinors projections for local layouts
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::LocalLayout > half_spinor_up_vector_t;
typedef Lattice::Lattice<Update::GaugeVector[2], Lattice::LocalLayout > half_spinor_down_vector_t;

typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > extended_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > standard_field_strength_lattice_t;
typedef Lattice::Lattice<Update::FermionicGroup[6], Lattice::LocalLayout > reduced_field_strength_lattice_t;

//Gauge transformation lattice equal for all layouts
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > extended_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > standard_matrix_lattice_t;
typedef Lattice::Lattice<Update::GaugeGroup, Lattice::LocalLayout > reduced_matrix_lattice_t;

//Gauge field configuration in the adjoint representations equal for all layouts
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::LocalLayout > extended_adjoint_lattice_t;
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::LocalLayout > standard_adjoint_lattice_t;
typedef Lattice::Lattice<Update::AdjointGroup[4], Lattice::LocalLayout > reduced_adjoint_lattice_t;

//Scalar field in the fundamental representation of the gauge group equal for all layouts
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::LocalLayout > extended_color_vector_t;
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::LocalLayout > standard_color_vector_t;
typedef Lattice::Lattice<Update::FundamentalVector, Lattice::LocalLayout > reduced_color_vector_t;

//Scalar field in the adjoint representation of the gauge group equal for all layouts
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::LocalLayout > extended_adjoint_complex_color_vector_t;
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::LocalLayout > standard_adjoint_complex_color_vector_t;
typedef Lattice::Lattice<Update::AdjointComplexVector, Lattice::LocalLayout > reduced_adjoint_complex_color_vector_t;

typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::LocalLayout > extended_adjoint_real_color_vector_t;
typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::LocalLayout > standard_adjoint_real_color_vector_t;
typedef Lattice::Lattice<Update::AdjointRealVector, Lattice::LocalLayout > reduced_adjoint_real_color_vector_t;

typedef Lattice::Lattice<int[4], Lattice::LocalLayout > extended_index_lattice_t;

#ifdef ALIGNED_OPT
typedef Lattice::ComplexVectorLattice<Update::real_t, Update::GaugeVector[4], Lattice::LocalLayout, 4, Update::diracVectorLength> reduced_soa_dirac_vector_t;
#endif

#endif

#ifdef ENABLE_MPI
typedef Lattice::Lattice<int, Lattice::MpiLayout<Lattice::ReducedStencil> > reduced_index_lattice_t;
#endif
#ifndef ENABLE_MPI
typedef Lattice::Lattice<int, Lattice::LocalLayout > reduced_index_lattice_t;
#endif

#endif

