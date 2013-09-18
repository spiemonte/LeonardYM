updater.o: updater.cpp
	$(CPP) $(CPPFLAGS) -c -o updater.o updater.cpp

LocalLayout.o: ./source/MPILattice/LocalLayout.h ./source/MPILattice/LocalLayout.cpp
	$(CPP) $(CPPFLAGS) -c -o LocalLayout.o ./source/MPILattice/LocalLayout.cpp

ReducedStencil.o: ./source/MPILattice/ReducedStencil.h ./source/MPILattice/ReducedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ReducedStencil.o ./source/MPILattice/ReducedStencil.cpp

StandardStencil.o: ./source/MPILattice/StandardStencil.h ./source/MPILattice/StandardStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o StandardStencil.o ./source/MPILattice/StandardStencil.cpp

ExtendedStencil.o: ./source/MPILattice/ExtendedStencil.h ./source/MPILattice/ExtendedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ExtendedStencil.o ./source/MPILattice/ExtendedStencil.cpp
	
AlgebraUtils.o: ./source/AlgebraUtils.h ./source/AlgebraUtils.cpp
	$(CPP) $(CPPFLAGS) -c -o AlgebraUtils.o ./source/AlgebraUtils.cpp

LatticeSweep.o: ./source/LatticeSweep.h ./source/LatticeSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o LatticeSweep.o ./source/LatticeSweep.cpp

Simulation.o: ./source/Simulation.h ./source/Simulation.cpp
	$(CPP) $(CPPFLAGS) -c -o Simulation.o ./source/Simulation.cpp

BiConjugateGradient.o: ./source/BiConjugateGradient.h ./source/BiConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o BiConjugateGradient.o ./source/BiConjugateGradient.cpp

ConjugateGradient.o: ./source/ConjugateGradient.h ./source/ConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o ConjugateGradient.o ./source/ConjugateGradient.cpp

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

SquareBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareBlockDiracWilsonOperator.h ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareBlockDiracWilsonOperator.o ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp

ComplementBlockDiracWilsonOperator.o: ./source/dirac_operators/ComplementBlockDiracWilsonOperator.h ./source/dirac_operators/ComplementBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ComplementBlockDiracWilsonOperator.o ./source/dirac_operators/ComplementBlockDiracWilsonOperator.cpp

SquareComplementBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.h ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareComplementBlockDiracWilsonOperator.o ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp

ImprovedDiracWilsonOperator.o: ./source/dirac_operators/ImprovedDiracWilsonOperator.h ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedDiracWilsonOperator.o ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp

SquareImprovedDiracWilsonOperator.o: ./source/dirac_operators/SquareImprovedDiracWilsonOperator.h ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o SquareImprovedDiracWilsonOperator.o ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp

MMMRMultishiftSolver.o: ./source/MMMRMultishiftSolver.h ./source/MMMRMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MMMRMultishiftSolver.o ./source/MMMRMultishiftSolver.cpp

MEMultishiftSolver.o: ./source/MEMultishiftSolver.h ./source/MEMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MEMultishiftSolver.o ./source/MEMultishiftSolver.cpp

MultishiftSolver.o: ./source/MultishiftSolver.h ./source/MultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o MultishiftSolver.o ./source/MultishiftSolver.cpp

ChronologicalMultishiftSolver.o: ./source/ChronologicalMultishiftSolver.h ./source/ChronologicalMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ChronologicalMultishiftSolver.o ./source/ChronologicalMultishiftSolver.cpp

DeflationInverter.o: ./source/DeflationInverter.h ./source/DeflationInverter.cpp
	$(CPP) $(CPPFLAGS) -c -o DeflationInverter.o ./source/DeflationInverter.cpp

WilsonGaugeAction.o: ./source/WilsonGaugeAction.h ./source/WilsonGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o WilsonGaugeAction.o ./source/WilsonGaugeAction.cpp

ImprovedGaugeAction.o: ./source/ImprovedGaugeAction.h ./source/ImprovedGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedGaugeAction.o ./source/ImprovedGaugeAction.cpp

GaugeForce.o: ./source/GaugeForce.h ./source/GaugeForce.cpp
	$(CPP) $(CPPFLAGS) -c -o GaugeForce.o ./source/GaugeForce.cpp

GaugeAction.o: ./source/GaugeAction.h ./source/GaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o GaugeAction.o ./source/GaugeAction.cpp

Integrate.o: ./source/Integrate.h ./source/Integrate.cpp
	$(CPP) $(CPPFLAGS) -c -o Integrate.o ./source/Integrate.cpp

LeapFrog.o: ./source/LeapFrog.h ./source/LeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o LeapFrog.o ./source/LeapFrog.cpp

FourthOrderLeapFrog.o: ./source/FourthOrderLeapFrog.h ./source/FourthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o FourthOrderLeapFrog.o ./source/FourthOrderLeapFrog.cpp

SixthOrderLeapFrog.o: ./source/SixthOrderLeapFrog.h ./source/SixthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o SixthOrderLeapFrog.o ./source/SixthOrderLeapFrog.cpp

OmelyanLeapFrog.o: ./source/OmelyanLeapFrog.h ./source/OmelyanLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o OmelyanLeapFrog.o ./source/OmelyanLeapFrog.cpp

Energy.o: ./source/Energy.h ./source/Energy.cpp
	$(CPP) $(CPPFLAGS) -c -o Energy.o ./source/Energy.cpp

Force.o: ./source/Force.h ./source/Force.cpp
	$(CPP) $(CPPFLAGS) -c -o Force.o ./source/Force.cpp

BandAction.o: ./source/BandAction.h ./source/BandAction.cpp
	$(CPP) $(CPPFLAGS) -c -o BandAction.o ./source/BandAction.cpp

BandTwoFlavorUpdater.o: ./source/BandTwoFlavorUpdater.h ./source/BandTwoFlavorUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o BandTwoFlavorUpdater.o ./source/BandTwoFlavorUpdater.cpp

HMCUpdater.o: ./source/HMCUpdater.h ./source/HMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o HMCUpdater.o ./source/HMCUpdater.cpp

FermionHMCUpdater.o: ./source/FermionHMCUpdater.h ./source/FermionHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionHMCUpdater.o ./source/FermionHMCUpdater.cpp

FermionicAction.o: ./source/FermionicAction.h ./source/FermionicAction.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionicAction.o ./source/FermionicAction.cpp

Polynomial.o: ./source/Polynomial.h ./source/Polynomial.cpp
	$(CPP) $(CPPFLAGS) -c -o Polynomial.o ./source/Polynomial.cpp

RationalApproximation.o: ./source/RationalApproximation.h ./source/RationalApproximation.cpp
	$(CPP) $(CPPFLAGS) -c -o RationalApproximation.o ./source/RationalApproximation.cpp

ReUnit.o: ./source/ReUnit.h ./source/ReUnit.cpp
	$(CPP) $(CPPFLAGS) -c -o ReUnit.o ./source/ReUnit.cpp

StoutSmearing.o: ./source/StoutSmearing.h ./source/StoutSmearing.cpp
	$(CPP) $(CPPFLAGS) -c -o StoutSmearing.o ./source/StoutSmearing.cpp

Gamma.o: ./source/Gamma.h ./source/Gamma.cpp
	$(CPP) $(CPPFLAGS) -c -o Gamma.o ./source/Gamma.cpp

PureGaugeUpdater.o: ./source/PureGaugeUpdater.h ./source/PureGaugeUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeUpdater.o ./source/PureGaugeUpdater.cpp

PureGaugeOverrelaxation.o: ./source/PureGaugeOverrelaxation.h ./source/PureGaugeOverrelaxation.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeOverrelaxation.o ./source/PureGaugeOverrelaxation.cpp

PureGaugeHMCUpdater.o: ./source/PureGaugeHMCUpdater.h ./source/PureGaugeHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeHMCUpdater.o ./source/PureGaugeHMCUpdater.cpp

Checkerboard.o: ./source/Checkerboard.h ./source/Checkerboard.cpp
	$(CPP) $(CPPFLAGS) -c -o Checkerboard.o ./source/Checkerboard.cpp

RandomSeed.o: ./source/RandomSeed.h ./source/RandomSeed.cpp
	$(CPP) $(CPPFLAGS) -c -o RandomSeed.o ./source/RandomSeed.cpp

Plaquette.o: ./source/Plaquette.h ./source/Plaquette.cpp
	$(CPP) $(CPPFLAGS) -c -o Plaquette.o ./source/Plaquette.cpp

PolyakovLoop.o: ./source/PolyakovLoop.h ./source/PolyakovLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o PolyakovLoop.o ./source/PolyakovLoop.cpp

WilsonLoop.o: ./source/WilsonLoop.h ./source/WilsonLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o WilsonLoop.o ./source/WilsonLoop.cpp

PureGaugeWilsonLoops.o: ./source/PureGaugeWilsonLoops.h ./source/PureGaugeWilsonLoops.cpp
	$(CPP) $(CPPFLAGS) -c -o PureGaugeWilsonLoops.o ./source/PureGaugeWilsonLoops.cpp 

StartGaugeConfiguration.o: ./source/StartGaugeConfiguration.h ./source/StartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o StartGaugeConfiguration.o ./source/StartGaugeConfiguration.cpp

ReadStartGaugeConfiguration.o: ./source/ReadStartGaugeConfiguration.h ./source/ReadStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ReadStartGaugeConfiguration.o ./source/ReadStartGaugeConfiguration.cpp

HotStartGaugeConfiguration.o: ./source/HotStartGaugeConfiguration.h ./source/HotStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o HotStartGaugeConfiguration.o ./source/HotStartGaugeConfiguration.cpp

ColdStartGaugeConfiguration.o: ./source/ColdStartGaugeConfiguration.h ./source/ColdStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ColdStartGaugeConfiguration.o ./source/ColdStartGaugeConfiguration.cpp

StochasticEstimator.o: ./source/StochasticEstimator.h ./source/StochasticEstimator.cpp
	$(CPP) $(CPPFLAGS) -c -o StochasticEstimator.o ./source/StochasticEstimator.cpp

MesonCorrelator.o: ./source/MesonCorrelator.h ./source/MesonCorrelator.cpp
	$(CPP) $(CPPFLAGS) -c -o MesonCorrelator.o ./source/MesonCorrelator.cpp

ChiralCondensate.o: ./source/ChiralCondensate.h ./source/ChiralCondensate.cpp
	$(CPP) $(CPPFLAGS) -c -o ChiralCondensate.o ./source/ChiralCondensate.cpp
	
GluinoGlue.o: ./source/GluinoGlue.h ./source/GluinoGlue.cpp
	$(CPP) $(CPPFLAGS) -c -o GluinoGlue.o ./source/GluinoGlue.cpp

FermionForce.o: ./source/FermionForce.h ./source/FermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o FermionForce.o ./source/FermionForce.cpp

DiracWilsonFermionForce.o: ./source/DiracWilsonFermionForce.h ./source/DiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracWilsonFermionForce.o ./source/DiracWilsonFermionForce.cpp

BlockDiracWilsonFermionForce.o: ./source/BlockDiracWilsonFermionForce.h ./source/BlockDiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o BlockDiracWilsonFermionForce.o ./source/BlockDiracWilsonFermionForce.cpp

ImprovedFermionForce.o: ./source/ImprovedFermionForce.h ./source/ImprovedFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ImprovedFermionForce.o ./source/ImprovedFermionForce.cpp

TestForce.o: ./source/TestForce.h ./source/TestForce.cpp
	$(CPP) $(CPPFLAGS) -c -o TestForce.o ./source/TestForce.cpp

TwoFlavorFermionAction.o: ./source/TwoFlavorFermionAction.h ./source/TwoFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorFermionAction.o ./source/TwoFlavorFermionAction.cpp

TwoFlavorQCDAction.o: ./source/TwoFlavorQCDAction.h ./source/TwoFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorQCDAction.o ./source/TwoFlavorQCDAction.cpp

TwoFlavorHMCUpdater.o: ./source/TwoFlavorHMCUpdater.h ./source/TwoFlavorHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o TwoFlavorHMCUpdater.o ./source/TwoFlavorHMCUpdater.cpp

NFlavorFermionAction.o: ./source/NFlavorFermionAction.h ./source/NFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorFermionAction.o ./source/NFlavorFermionAction.cpp

NFlavorQCDAction.o: ./source/NFlavorQCDAction.h ./source/NFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorQCDAction.o ./source/NFlavorQCDAction.cpp

NFlavorQCDUpdater.o : ./source/NFlavorQCDUpdater.h ./source/NFlavorQCDUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorQCDUpdater.o ./source/NFlavorQCDUpdater.cpp

NFlavorBlockUpdater.o : ./source/NFlavorBlockUpdater.h ./source/NFlavorBlockUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o NFlavorBlockUpdater.o ./source/NFlavorBlockUpdater.cpp

DiracEigenSolver.o: ./source/DiracEigenSolver.h ./source/DiracEigenSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o DiracEigenSolver.o ./source/DiracEigenSolver.cpp

Eigenvalues.o: ./source/Eigenvalues.h ./source/Eigenvalues.cpp
	$(CPP) $(CPPFLAGS) -c -o Eigenvalues.o ./source/Eigenvalues.cpp

GlobalOutput.o: ./source/GlobalOutput.h ./source/GlobalOutput.cpp
	$(CPP) $(CPPFLAGS) -c -o GlobalOutput.o ./source/GlobalOutput.cpp

OutputSweep.o: ./source/OutputSweep.h ./source/OutputSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o OutputSweep.o ./source/OutputSweep.cpp

Glueball.o: ./source/Glueball.h ./source/Glueball.cpp
	$(CPP) $(CPPFLAGS) -c -o Glueball.o ./source/Glueball.cpp

TestCommunication.o: ./source/TestCommunication.h ./source/TestCommunication.cpp
	$(CPP) $(CPPFLAGS) -c -o TestCommunication.o ./source/TestCommunication.cpp

TestLinearAlgebra.o: ./source/TestLinearAlgebra.h ./source/TestLinearAlgebra.cpp
	$(CPP) $(CPPFLAGS) -c -o TestLinearAlgebra.o ./source/TestLinearAlgebra.cpp

StorageParameters.o: ./source/StorageParameters.h ./source/StorageParameters.cpp
	$(CPP) $(CPPFLAGS) -c -o StorageParameters.o ./source/StorageParameters.cpp
