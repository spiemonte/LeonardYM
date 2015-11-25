updater.o: updater.cpp
	$(CPP) $(CPPFLAGS) -c -o updater.o updater.cpp

Environment.o: ./source/Environment.cpp ./source/Environment.h
	$(CPP) $(CPPFLAGS) -c -o Environment.o ./source/Environment.cpp

LocalLayout.o: ./source/MPILattice/LocalLayout.h ./source/MPILattice/LocalLayout.cpp
	$(CPP) $(CPPFLAGS) -c -o LocalLayout.o ./source/MPILattice/LocalLayout.cpp

ReducedStencil.o: ./source/MPILattice/ReducedStencil.h ./source/MPILattice/ReducedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ReducedStencil.o ./source/MPILattice/ReducedStencil.cpp

StandardStencil.o: ./source/MPILattice/StandardStencil.h ./source/MPILattice/StandardStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o StandardStencil.o ./source/MPILattice/StandardStencil.cpp

ExtendedStencil.o: ./source/MPILattice/ExtendedStencil.h ./source/MPILattice/ExtendedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ExtendedStencil.o ./source/MPILattice/ExtendedStencil.cpp
	
AlgebraUtils.o: ./source/algebra_utils/AlgebraUtils.h ./source/algebra_utils/AlgebraUtils.cpp
	$(CPP) $(CPPFLAGS) -c -o AlgebraUtils.o ./source/algebra_utils/AlgebraUtils.cpp

LatticeSweep.o: ./source/LatticeSweep.h ./source/LatticeSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o LatticeSweep.o ./source/LatticeSweep.cpp

Simulation.o: ./source/Simulation.h ./source/Simulation.cpp
	$(CPP) $(CPPFLAGS) -c -o Simulation.o ./source/Simulation.cpp

BiConjugateGradient.o: ./source/inverters/BiConjugateGradient.h ./source/inverters/BiConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o BiConjugateGradient.o ./source/inverters/BiConjugateGradient.cpp

GMRESR.o: ./source/inverters/GMRESR.h ./source/inverters/GMRESR.cpp
	$(CPP) $(CPPFLAGS) -c -o GMRESR.o ./source/inverters/GMRESR.cpp

ConjugateGradient.o: ./source/inverters/ConjugateGradient.h ./source/inverters/ConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o ConjugateGradient.o ./source/inverters/ConjugateGradient.cpp

DiracOperator.o: ./source/dirac_operators/DiracOperator.h ./source/dirac_operators/DiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracOperator.o ./source/dirac_operators/DiracOperator.cpp

BasicDiracWilsonOperator.o: ./source/dirac_operators/BasicDiracWilsonOperator.h ./source/dirac_operators/BasicDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o BasicDiracWilsonOperator.o ./source/dirac_operators/BasicDiracWilsonOperator.cpp

BasicSquareDiracWilsonOperator.o: ./source/dirac_operators/BasicSquareDiracWilsonOperator.h ./source/dirac_operators/BasicSquareDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o BasicSquareDiracWilsonOperator.o ./source/dirac_operators/BasicSquareDiracWilsonOperator.cpp

DiracWilsonOperator.o: ./source/dirac_operators/DiracWilsonOperator.h ./source/dirac_operators/DiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracWilsonOperator.o ./source/dirac_operators/DiracWilsonOperator.cpp

SquareDiracWilsonOperator.o: ./source/dirac_operators/SquareDiracWilsonOperator.h ./source/dirac_operators/SquareDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareDiracWilsonOperator.o ./source/dirac_operators/SquareDiracWilsonOperator.cpp

BlockDiracOperator.o: ./source/dirac_operators/BlockDiracOperator.h ./source/dirac_operators/BlockDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o BlockDiracOperator.o ./source/dirac_operators/BlockDiracOperator.cpp

BlockDiracWilsonOperator.o: ./source/dirac_operators/BlockDiracWilsonOperator.h ./source/dirac_operators/BlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o BlockDiracWilsonOperator.o ./source/dirac_operators/BlockDiracWilsonOperator.cpp

BlockImprovedDiracWilsonOperator.o: ./source/dirac_operators/BlockImprovedDiracWilsonOperator.h ./source/dirac_operators/BlockImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o BlockImprovedDiracWilsonOperator.o ./source/dirac_operators/BlockImprovedDiracWilsonOperator.cpp

SquareBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareBlockDiracWilsonOperator.h ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareBlockDiracWilsonOperator.o ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp

SquareTwistedDiracOperator.o: ./source/dirac_operators/SquareTwistedDiracOperator.h ./source/dirac_operators/SquareTwistedDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareTwistedDiracOperator.o ./source/dirac_operators/SquareTwistedDiracOperator.cpp

TwistedDiracOperator.o: ./source/dirac_operators/TwistedDiracOperator.h ./source/dirac_operators/TwistedDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o TwistedDiracOperator.o ./source/dirac_operators/TwistedDiracOperator.cpp

ComplementBlockDiracWilsonOperator.o: ./source/dirac_operators/ComplementBlockDiracWilsonOperator.h ./source/dirac_operators/ComplementBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ComplementBlockDiracWilsonOperator.o ./source/dirac_operators/ComplementBlockDiracWilsonOperator.cpp

SquareComplementBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.h ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareComplementBlockDiracWilsonOperator.o ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp

SquareComplementBlockDiracOperator.o: ./source/dirac_operators/SquareComplementBlockDiracOperator.h ./source/dirac_operators/SquareComplementBlockDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareComplementBlockDiracOperator.o ./source/dirac_operators/SquareComplementBlockDiracOperator.cpp

ImprovedDiracWilsonOperator.o: ./source/dirac_operators/ImprovedDiracWilsonOperator.h ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedDiracWilsonOperator.o ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp

SquareImprovedDiracWilsonOperator.o: ./source/dirac_operators/SquareImprovedDiracWilsonOperator.h ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareImprovedDiracWilsonOperator.o ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp

MMMRMultishiftSolver.o: ./source/inverters/MMMRMultishiftSolver.h ./source/inverters/MMMRMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MMMRMultishiftSolver.o ./source/inverters/MMMRMultishiftSolver.cpp

MEMultishiftSolver.o: ./source/inverters/MEMultishiftSolver.h ./source/inverters/MEMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MEMultishiftSolver.o ./source/inverters/MEMultishiftSolver.cpp

MultishiftSolver.o: ./source/inverters/MultishiftSolver.h ./source/inverters/MultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MultishiftSolver.o ./source/inverters/MultishiftSolver.cpp

ChronologicalMultishiftSolver.o: ./source/inverters/ChronologicalMultishiftSolver.h ./source/inverters/ChronologicalMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ChronologicalMultishiftSolver.o ./source/inverters/ChronologicalMultishiftSolver.cpp

DeflationInverter.o: ./source/inverters/DeflationInverter.h ./source/inverters/DeflationInverter.cpp
	$(CPP) $(CPPFLAGS) -c -o DeflationInverter.o ./source/inverters/DeflationInverter.cpp

WilsonGaugeAction.o: ./source/actions/WilsonGaugeAction.h ./source/actions/WilsonGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o WilsonGaugeAction.o ./source/actions/WilsonGaugeAction.cpp

ImprovedGaugeAction.o: ./source/actions/ImprovedGaugeAction.h ./source/actions/ImprovedGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedGaugeAction.o ./source/actions/ImprovedGaugeAction.cpp

GaugeForce.o: ./source/hmc_forces/GaugeForce.h ./source/hmc_forces/GaugeForce.cpp
	$(CPP) $(CPPFLAGS) -c -o GaugeForce.o ./source/hmc_forces/GaugeForce.cpp

GaugeAction.o: ./source/actions/GaugeAction.h ./source/actions/GaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o GaugeAction.o ./source/actions/GaugeAction.cpp

Integrate.o: ./source/hmc_integrators/Integrate.h ./source/hmc_integrators/Integrate.cpp
	$(CPP) $(CPPFLAGS) -c -o Integrate.o ./source/hmc_integrators/Integrate.cpp

LeapFrog.o: ./source/hmc_integrators/LeapFrog.h ./source/hmc_integrators/LeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o LeapFrog.o ./source/hmc_integrators/LeapFrog.cpp

FourthOrderLeapFrog.o: ./source/hmc_integrators/FourthOrderLeapFrog.h ./source/hmc_integrators/FourthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o FourthOrderLeapFrog.o ./source/hmc_integrators/FourthOrderLeapFrog.cpp

SixthOrderLeapFrog.o: ./source/hmc_integrators/SixthOrderLeapFrog.h ./source/hmc_integrators/SixthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o SixthOrderLeapFrog.o ./source/hmc_integrators/SixthOrderLeapFrog.cpp

OmelyanLeapFrog.o: ./source/hmc_integrators/OmelyanLeapFrog.h ./source/hmc_integrators/OmelyanLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o OmelyanLeapFrog.o ./source/hmc_integrators/OmelyanLeapFrog.cpp

FourthOmelyanLeapFrog.o: ./source/hmc_integrators/FourthOmelyanLeapFrog.h ./source/hmc_integrators/FourthOmelyanLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o FourthOmelyanLeapFrog.o ./source/hmc_integrators/FourthOmelyanLeapFrog.cpp

Energy.o: ./source/actions/Energy.h ./source/actions/Energy.cpp
	$(CPP) $(CPPFLAGS) -c -o Energy.o ./source/actions/Energy.cpp

Force.o: ./source/hmc_forces/Force.h ./source/hmc_forces/Force.cpp
	$(CPP) $(CPPFLAGS) -c -o Force.o ./source/hmc_forces/Force.cpp

BandAction.o: ./source/actions/BandAction.h ./source/actions/BandAction.cpp
	$(CPP) $(CPPFLAGS) -c -o BandAction.o ./source/actions/BandAction.cpp

BandTwoFlavorUpdater.o: ./source/hmc_updaters/BandTwoFlavorUpdater.h ./source/hmc_updaters/BandTwoFlavorUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o BandTwoFlavorUpdater.o ./source/hmc_updaters/BandTwoFlavorUpdater.cpp

HMCUpdater.o: ./source/hmc_updaters/HMCUpdater.h ./source/hmc_updaters/HMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o HMCUpdater.o ./source/hmc_updaters/HMCUpdater.cpp

FermionHMCUpdater.o: ./source/hmc_updaters/FermionHMCUpdater.h ./source/hmc_updaters/FermionHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionHMCUpdater.o ./source/hmc_updaters/FermionHMCUpdater.cpp

FermionicAction.o: ./source/actions/FermionicAction.h ./source/actions/FermionicAction.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionicAction.o ./source/actions/FermionicAction.cpp

Polynomial.o: ./source/dirac_functions/Polynomial.h ./source/dirac_functions/Polynomial.cpp
	$(CPP) $(CPPFLAGS) -c -o Polynomial.o ./source/dirac_functions/Polynomial.cpp

RationalApproximation.o: ./source/dirac_functions/RationalApproximation.h ./source/dirac_functions/RationalApproximation.cpp
	$(CPP) $(CPPFLAGS) -c -o RationalApproximation.o ./source/dirac_functions/RationalApproximation.cpp

ReUnit.o: ./source/utils/ReUnit.h ./source/utils/ReUnit.cpp
	$(CPP) $(CPPFLAGS) -c -o ReUnit.o ./source/utils/ReUnit.cpp

StoutSmearing.o: ./source/utils/StoutSmearing.h ./source/utils/StoutSmearing.cpp
	$(CPP) $(CPPFLAGS) -c -o StoutSmearing.o ./source/utils/StoutSmearing.cpp

Gamma.o: ./source/utils/Gamma.h ./source/utils/Gamma.cpp
	$(CPP) $(CPPFLAGS) -c -o Gamma.o ./source/utils/Gamma.cpp

PureGaugeUpdater.o: ./source/pure_gauge/PureGaugeUpdater.h ./source/pure_gauge/PureGaugeUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeUpdater.o ./source/pure_gauge/PureGaugeUpdater.cpp

PureGaugeOverrelaxation.o: ./source/pure_gauge/PureGaugeOverrelaxation.h ./source/pure_gauge/PureGaugeOverrelaxation.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeOverrelaxation.o ./source/pure_gauge/PureGaugeOverrelaxation.cpp

PureGaugeHMCUpdater.o: ./source/hmc_updaters/PureGaugeHMCUpdater.h ./source/hmc_updaters/PureGaugeHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeHMCUpdater.o ./source/hmc_updaters/PureGaugeHMCUpdater.cpp

Checkerboard.o: ./source/pure_gauge/Checkerboard.h ./source/pure_gauge/Checkerboard.cpp
	$(CPP) $(CPPFLAGS) -c -o Checkerboard.o ./source/pure_gauge/Checkerboard.cpp

RandomSeed.o: ./source/utils/RandomSeed.h ./source/utils/RandomSeed.cpp
	$(CPP) $(CPPFLAGS) -c -o RandomSeed.o ./source/utils/RandomSeed.cpp

Plaquette.o: ./source/wilson_loops/Plaquette.h ./source/wilson_loops/Plaquette.cpp
	$(CPP) $(CPPFLAGS) -c -o Plaquette.o ./source/wilson_loops/Plaquette.cpp

PolyakovLoop.o: ./source/polyakov_loops/PolyakovLoop.h ./source/polyakov_loops/PolyakovLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o PolyakovLoop.o ./source/polyakov_loops/PolyakovLoop.cpp

PolyakovLoopEigenvalues.o: ./source/polyakov_loops/PolyakovLoopEigenvalues.h ./source/polyakov_loops/PolyakovLoopEigenvalues.cpp
	$(CPP) $(CPPFLAGS) -c -o PolyakovLoopEigenvalues.o ./source/polyakov_loops/PolyakovLoopEigenvalues.cpp

PolyakovLoopCorrelator.o: ./source/polyakov_loops/PolyakovLoopCorrelator.h ./source/polyakov_loops/PolyakovLoopCorrelator.cpp
	$(CPP) $(CPPFLAGS) -c -o PolyakovLoopCorrelator.o ./source/polyakov_loops/PolyakovLoopCorrelator.cpp

AdjointPolyakovLoop.o: ./source/polyakov_loops/AdjointPolyakovLoop.h ./source/polyakov_loops/AdjointPolyakovLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o AdjointPolyakovLoop.o ./source/polyakov_loops/AdjointPolyakovLoop.cpp

WilsonLoop.o: ./source/wilson_loops/WilsonLoop.h ./source/wilson_loops/WilsonLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o WilsonLoop.o ./source/wilson_loops/WilsonLoop.cpp

PureGaugeWilsonLoops.o: ./source/pure_gauge/PureGaugeWilsonLoops.h ./source/pure_gauge/PureGaugeWilsonLoops.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeWilsonLoops.o ./source/pure_gauge/PureGaugeWilsonLoops.cpp 

StartGaugeConfiguration.o: ./source/starters/StartGaugeConfiguration.h ./source/starters/StartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o StartGaugeConfiguration.o ./source/starters/StartGaugeConfiguration.cpp

ReadStartGaugeConfiguration.o: ./source/starters/ReadStartGaugeConfiguration.h ./source/starters/ReadStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ReadStartGaugeConfiguration.o ./source/starters/ReadStartGaugeConfiguration.cpp

HotStartGaugeConfiguration.o: ./source/starters/HotStartGaugeConfiguration.h ./source/starters/HotStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o HotStartGaugeConfiguration.o ./source/starters/HotStartGaugeConfiguration.cpp

ColdStartGaugeConfiguration.o: ./source/starters/ColdStartGaugeConfiguration.h ./source/starters/ColdStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ColdStartGaugeConfiguration.o ./source/starters/ColdStartGaugeConfiguration.cpp

StochasticEstimator.o: ./source/fermion_measurements/StochasticEstimator.h ./source/fermion_measurements/StochasticEstimator.cpp
	$(CPP) $(CPPFLAGS) -c -o StochasticEstimator.o ./source/fermion_measurements/StochasticEstimator.cpp

MesonCorrelator.o: ./source/correlators/MesonCorrelator.h ./source/correlators/MesonCorrelator.cpp
	$(CPP) $(CPPFLAGS) -c -o MesonCorrelator.o ./source/correlators/MesonCorrelator.cpp

ChiralCondensate.o: ./source/fermion_measurements/ChiralCondensate.h ./source/fermion_measurements/ChiralCondensate.cpp
	$(CPP) $(CPPFLAGS) -c -o ChiralCondensate.o ./source/fermion_measurements/ChiralCondensate.cpp

NPRVertex.o: ./source/fermion_measurements/NPRVertex.h ./source/fermion_measurements/NPRVertex.cpp
	$(CPP) $(CPPFLAGS) -c -o NPRVertex.o ./source/fermion_measurements/NPRVertex.cpp
	
GluinoGlue.o: ./source/correlators/GluinoGlue.h ./source/correlators/GluinoGlue.cpp
	$(CPP) $(CPPFLAGS) -c -o GluinoGlue.o ./source/correlators/GluinoGlue.cpp

FermionForce.o: ./source/hmc_forces/FermionForce.h ./source/hmc_forces/FermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionForce.o ./source/hmc_forces/FermionForce.cpp

SmearingForce.o: ./source/hmc_forces/SmearingForce.h ./source/hmc_forces/SmearingForce.cpp
	$(CPP) $(CPPFLAGS) -c -o SmearingForce.o ./source/hmc_forces/SmearingForce.cpp

DiracWilsonFermionForce.o: ./source/hmc_forces/DiracWilsonFermionForce.h ./source/hmc_forces/DiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracWilsonFermionForce.o ./source/hmc_forces/DiracWilsonFermionForce.cpp

BlockDiracWilsonFermionForce.o: ./source/hmc_forces/BlockDiracWilsonFermionForce.h ./source/hmc_forces/BlockDiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o BlockDiracWilsonFermionForce.o ./source/hmc_forces/BlockDiracWilsonFermionForce.cpp

ImprovedFermionForce.o: ./source/hmc_forces/ImprovedFermionForce.h ./source/hmc_forces/ImprovedFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedFermionForce.o ./source/hmc_forces/ImprovedFermionForce.cpp

TestForce.o: ./source/hmc_forces/TestForce.h ./source/hmc_forces/TestForce.cpp
	$(CPP) $(CPPFLAGS) -c -o TestForce.o ./source/hmc_forces/TestForce.cpp

TwoFlavorFermionAction.o: ./source/actions/TwoFlavorFermionAction.h ./source/actions/TwoFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorFermionAction.o ./source/actions/TwoFlavorFermionAction.cpp

TwoFlavorQCDAction.o: ./source/actions/TwoFlavorQCDAction.h ./source/actions/TwoFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorQCDAction.o ./source/actions/TwoFlavorQCDAction.cpp

TwoFlavorHMCUpdater.o: ./source/hmc_updaters/TwoFlavorHMCUpdater.h ./source/hmc_updaters/TwoFlavorHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorHMCUpdater.o ./source/hmc_updaters/TwoFlavorHMCUpdater.cpp

NFlavorFermionAction.o: ./source/actions/NFlavorFermionAction.h ./source/actions/NFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorFermionAction.o ./source/actions/NFlavorFermionAction.cpp

NFlavorQCDAction.o: ./source/actions/NFlavorQCDAction.h ./source/actions/NFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorQCDAction.o ./source/actions/NFlavorQCDAction.cpp

NFlavorQCDUpdater.o : ./source/hmc_updaters/NFlavorQCDUpdater.h ./source/hmc_updaters/NFlavorQCDUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorQCDUpdater.o ./source/hmc_updaters/NFlavorQCDUpdater.cpp

MultiStepNFlavorQCDUpdater.o : ./source/hmc_updaters/MultiStepNFlavorQCDUpdater.h ./source/hmc_updaters/MultiStepNFlavorQCDUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o MultiStepNFlavorQCDUpdater.o ./source/hmc_updaters/MultiStepNFlavorQCDUpdater.cpp

TwistedMultiStepNFlavorQCDUpdater.o : ./source/hmc_updaters/TwistedMultiStepNFlavorQCDUpdater.h ./source/hmc_updaters/TwistedMultiStepNFlavorQCDUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o TwistedMultiStepNFlavorQCDUpdater.o ./source/hmc_updaters/TwistedMultiStepNFlavorQCDUpdater.cpp

NFlavorBlockUpdater.o : ./source/hmc_updaters/NFlavorBlockUpdater.h ./source/hmc_updaters/NFlavorBlockUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorBlockUpdater.o ./source/hmc_updaters/NFlavorBlockUpdater.cpp

DiracEigenSolver.o: ./source/fermion_measurements/DiracEigenSolver.h ./source/fermion_measurements/DiracEigenSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracEigenSolver.o ./source/fermion_measurements/DiracEigenSolver.cpp

SingletOperators.o: ./source/fermion_measurements/SingletOperators.h ./source/fermion_measurements/SingletOperators.cpp
	$(CPP) $(CPPFLAGS) -c -o SingletOperators.o ./source/fermion_measurements/SingletOperators.cpp

Eigenvalues.o: ./source/fermion_measurements/Eigenvalues.h ./source/fermion_measurements/Eigenvalues.cpp
	$(CPP) $(CPPFLAGS) -c -o Eigenvalues.o ./source/fermion_measurements/Eigenvalues.cpp

GlobalOutput.o: ./source/io/GlobalOutput.h ./source/io/GlobalOutput.cpp
	$(CPP) $(CPPFLAGS) -c -o GlobalOutput.o ./source/io/GlobalOutput.cpp

OutputSweep.o: ./source/io/OutputSweep.h ./source/io/OutputSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o OutputSweep.o ./source/io/OutputSweep.cpp

Glueball.o: ./source/correlators/Glueball.h ./source/correlators/Glueball.cpp
	$(CPP) $(CPPFLAGS) -c -o Glueball.o ./source/correlators/Glueball.cpp

TestCommunication.o: ./source/tests/TestCommunication.h ./source/tests/TestCommunication.cpp
	$(CPP) $(CPPFLAGS) -c -o TestCommunication.o ./source/tests/TestCommunication.cpp

TestLinearAlgebra.o: ./source/tests/TestLinearAlgebra.h ./source/tests/TestLinearAlgebra.cpp
	$(CPP) $(CPPFLAGS) -c -o TestLinearAlgebra.o ./source/tests/TestLinearAlgebra.cpp

TestSpeedDiracOperators.o: ./source/tests/TestSpeedDiracOperators.h ./source/tests/TestSpeedDiracOperators.cpp
	$(CPP) $(CPPFLAGS) -c -o TestSpeedDiracOperators.o ./source/tests/TestSpeedDiracOperators.cpp

StorageParameters.o: ./source/io/StorageParameters.h ./source/io/StorageParameters.cpp
	$(CPP) $(CPPFLAGS) -c -o StorageParameters.o ./source/io/StorageParameters.cpp

ReadGaugeConfiguration.o: ./source/io/ReadGaugeConfiguration.h ./source/io/ReadGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ReadGaugeConfiguration.o ./source/io/ReadGaugeConfiguration.cpp
	
WilsonFlow.o: ./source/wilson_flow/WilsonFlow.h ./source/wilson_flow/WilsonFlow.cpp
	$(CPP) $(CPPFLAGS) -c -o WilsonFlow.o ./source/wilson_flow/WilsonFlow.cpp

GaugeEnergy.o: ./source/actions/GaugeEnergy.h ./source/actions/GaugeEnergy.cpp
	$(CPP) $(CPPFLAGS) -c -o GaugeEnergy.o ./source/actions/GaugeEnergy.cpp

MultiGridDiracOperator.o: ./source/multigrid_solver/MultiGridDiracOperator.h ./source/multigrid_solver/MultiGridDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o MultiGridDiracOperator.o ./source/multigrid_solver/MultiGridDiracOperator.cpp

MultiGridSolver.o: ./source/multigrid_solver/MultiGridSolver.h ./source/multigrid_solver/MultiGridSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MultiGridSolver.o ./source/multigrid_solver/MultiGridSolver.cpp

SAPPreconditioner.o: ./source/dirac_operators/SAPPreconditioner.h ./source/dirac_operators/SAPPreconditioner.cpp
	$(CPP) $(CPPFLAGS) -c -o SAPPreconditioner.o ./source/dirac_operators/SAPPreconditioner.cpp
	
MultiGridOperator.o: ./source/dirac_operators/MultiGridOperator.h ./source/dirac_operators/MultiGridOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o MultiGridOperator.o ./source/dirac_operators/MultiGridOperator.cpp

