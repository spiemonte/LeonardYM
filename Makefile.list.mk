./build/main.o: ./source/main.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/main.o ./source/main.cpp

./build/Environment.o: ./source/Environment.cpp ./source/Environment.h
	$(CPP) $(CPPFLAGS) -c -o ./build/Environment.o ./source/Environment.cpp

./build/LocalLayout.o: ./source/MPILattice/LocalLayout.h ./source/MPILattice/LocalLayout.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LocalLayout.o ./source/MPILattice/LocalLayout.cpp

./build/ReducedStencil.o: ./source/MPILattice/ReducedStencil.h ./source/MPILattice/ReducedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ReducedStencil.o ./source/MPILattice/ReducedStencil.cpp

./build/StandardStencil.o: ./source/MPILattice/StandardStencil.h ./source/MPILattice/StandardStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/StandardStencil.o ./source/MPILattice/StandardStencil.cpp

./build/ExtendedStencil.o: ./source/MPILattice/ExtendedStencil.h ./source/MPILattice/ExtendedStencil.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ExtendedStencil.o ./source/MPILattice/ExtendedStencil.cpp
	
./build/AlgebraUtils.o: ./source/algebra_utils/AlgebraUtils.h ./source/algebra_utils/AlgebraUtils.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/AlgebraUtils.o ./source/algebra_utils/AlgebraUtils.cpp

./build/LatticeSweep.o: ./source/LatticeSweep.h ./source/LatticeSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LatticeSweep.o ./source/LatticeSweep.cpp

./build/Simulation.o: ./source/Simulation.h ./source/Simulation.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Simulation.o ./source/Simulation.cpp

./build/BiConjugateGradient.o: ./source/inverters/BiConjugateGradient.h ./source/inverters/BiConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BiConjugateGradient.o ./source/inverters/BiConjugateGradient.cpp

./build/PreconditionedBiCGStab.o: ./source/inverters/PreconditionedBiCGStab.h ./source/inverters/PreconditionedBiCGStab.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PreconditionedBiCGStab.o ./source/inverters/PreconditionedBiCGStab.cpp

./build/GMRESR.o: ./source/inverters/GMRESR.h ./source/inverters/GMRESR.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GMRESR.o ./source/inverters/GMRESR.cpp

./build/ConjugateGradient.o: ./source/inverters/ConjugateGradient.h ./source/inverters/ConjugateGradient.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ConjugateGradient.o ./source/inverters/ConjugateGradient.cpp

./build/DiracOperator.o: ./source/dirac_operators/DiracOperator.h ./source/dirac_operators/DiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/DiracOperator.o ./source/dirac_operators/DiracOperator.cpp

./build/Propagator.o: ./source/dirac_operators/Propagator.h ./source/dirac_operators/Propagator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Propagator.o ./source/dirac_operators/Propagator.cpp

./build/BasicDiracWilsonOperator.o: ./source/dirac_operators/BasicDiracWilsonOperator.h ./source/dirac_operators/BasicDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BasicDiracWilsonOperator.o ./source/dirac_operators/BasicDiracWilsonOperator.cpp

./build/BasicSquareDiracWilsonOperator.o: ./source/dirac_operators/BasicSquareDiracWilsonOperator.h ./source/dirac_operators/BasicSquareDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BasicSquareDiracWilsonOperator.o ./source/dirac_operators/BasicSquareDiracWilsonOperator.cpp

./build/DiracWilsonOperator.o: ./source/dirac_operators/DiracWilsonOperator.h ./source/dirac_operators/DiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/DiracWilsonOperator.o ./source/dirac_operators/DiracWilsonOperator.cpp

./build/SquareDiracWilsonOperator.o: ./source/dirac_operators/SquareDiracWilsonOperator.h ./source/dirac_operators/SquareDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareDiracWilsonOperator.o ./source/dirac_operators/SquareDiracWilsonOperator.cpp

./build/OverlapOperator.o: ./source/dirac_operators/OverlapOperator.h ./source/dirac_operators/OverlapOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/OverlapOperator.o ./source/dirac_operators/OverlapOperator.cpp

./build/SquareOverlapOperator.o: ./source/dirac_operators/SquareOverlapOperator.h ./source/dirac_operators/SquareOverlapOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareOverlapOperator.o ./source/dirac_operators/SquareOverlapOperator.cpp

./build/ExactOverlapOperator.o: ./source/dirac_operators/ExactOverlapOperator.h ./source/dirac_operators/ExactOverlapOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ExactOverlapOperator.o ./source/dirac_operators/ExactOverlapOperator.cpp

./build/BlockDiracOperator.o: ./source/dirac_operators/BlockDiracOperator.h ./source/dirac_operators/BlockDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BlockDiracOperator.o ./source/dirac_operators/BlockDiracOperator.cpp

./build/BlockDiracWilsonOperator.o: ./source/dirac_operators/BlockDiracWilsonOperator.h ./source/dirac_operators/BlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BlockDiracWilsonOperator.o ./source/dirac_operators/BlockDiracWilsonOperator.cpp

./build/BlockImprovedDiracWilsonOperator.o: ./source/dirac_operators/BlockImprovedDiracWilsonOperator.h ./source/dirac_operators/BlockImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BlockImprovedDiracWilsonOperator.o ./source/dirac_operators/BlockImprovedDiracWilsonOperator.cpp

./build/SquareBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareBlockDiracWilsonOperator.h ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareBlockDiracWilsonOperator.o ./source/dirac_operators/SquareBlockDiracWilsonOperator.cpp

./build/SquareTwistedDiracOperator.o: ./source/dirac_operators/SquareTwistedDiracOperator.h ./source/dirac_operators/SquareTwistedDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareTwistedDiracOperator.o ./source/dirac_operators/SquareTwistedDiracOperator.cpp

./build/TwistedDiracOperator.o: ./source/dirac_operators/TwistedDiracOperator.h ./source/dirac_operators/TwistedDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TwistedDiracOperator.o ./source/dirac_operators/TwistedDiracOperator.cpp

./build/ComplementBlockDiracOperator.o: ./source/dirac_operators/ComplementBlockDiracOperator.h ./source/dirac_operators/ComplementBlockDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ComplementBlockDiracOperator.o ./source/dirac_operators/ComplementBlockDiracOperator.cpp

./build/SquareComplementBlockDiracWilsonOperator.o: ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.h ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareComplementBlockDiracWilsonOperator.o ./source/dirac_operators/SquareComplementBlockDiracWilsonOperator.cpp

./build/HoppingOperator.o: ./source/dirac_operators/HoppingOperator.h ./source/dirac_operators/HoppingOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/HoppingOperator.o ./source/dirac_operators/HoppingOperator.cpp

./build/GammaOperators.o: ./source/dirac_operators/GammaOperators.h ./source/dirac_operators/GammaOperators.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GammaOperators.o ./source/dirac_operators/GammaOperators.cpp

./build/SquareComplementBlockDiracOperator.o: ./source/dirac_operators/SquareComplementBlockDiracOperator.h ./source/dirac_operators/SquareComplementBlockDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareComplementBlockDiracOperator.o ./source/dirac_operators/SquareComplementBlockDiracOperator.cpp

./build/ImprovedDiracWilsonOperator.o: ./source/dirac_operators/ImprovedDiracWilsonOperator.h ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ImprovedDiracWilsonOperator.o ./source/dirac_operators/ImprovedDiracWilsonOperator.cpp

./build/EvenOddImprovedDiracWilsonOperator.o: ./source/dirac_operators/EvenOddImprovedDiracWilsonOperator.h ./source/dirac_operators/EvenOddImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/EvenOddImprovedDiracWilsonOperator.o ./source/dirac_operators/EvenOddImprovedDiracWilsonOperator.cpp

./build/SquareEvenOddImprovedDiracWilsonOperator.o: ./source/dirac_operators/SquareEvenOddImprovedDiracWilsonOperator.h ./source/dirac_operators/EvenOddImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareEvenOddImprovedDiracWilsonOperator.o ./source/dirac_operators/SquareEvenOddImprovedDiracWilsonOperator.cpp

./build/SquareImprovedDiracWilsonOperator.o: ./source/dirac_operators/SquareImprovedDiracWilsonOperator.h ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SquareImprovedDiracWilsonOperator.o ./source/dirac_operators/SquareImprovedDiracWilsonOperator.cpp

./build/MMMRMultishiftSolver.o: ./source/inverters/MMMRMultishiftSolver.h ./source/inverters/MMMRMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MMMRMultishiftSolver.o ./source/inverters/MMMRMultishiftSolver.cpp

./build/MEMultishiftSolver.o: ./source/inverters/MEMultishiftSolver.h ./source/inverters/MEMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MEMultishiftSolver.o ./source/inverters/MEMultishiftSolver.cpp

./build/MultishiftSolver.o: ./source/inverters/MultishiftSolver.h ./source/inverters/MultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MultishiftSolver.o ./source/inverters/MultishiftSolver.cpp

./build/ChronologicalMultishiftSolver.o: ./source/inverters/ChronologicalMultishiftSolver.h ./source/inverters/ChronologicalMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ChronologicalMultishiftSolver.o ./source/inverters/ChronologicalMultishiftSolver.cpp

./build/MultiGridMEMultishiftSolver.o: ./source/inverters/MultiGridMEMultishiftSolver.h ./source/inverters/MultiGridMEMultishiftSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridMEMultishiftSolver.o ./source/inverters/MultiGridMEMultishiftSolver.cpp

./build/DeflationInverter.o: ./source/inverters/DeflationInverter.h ./source/inverters/DeflationInverter.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/DeflationInverter.o ./source/inverters/DeflationInverter.cpp

./build/WilsonGaugeAction.o: ./source/actions/WilsonGaugeAction.h ./source/actions/WilsonGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/WilsonGaugeAction.o ./source/actions/WilsonGaugeAction.cpp

./build/ImprovedGaugeAction.o: ./source/actions/ImprovedGaugeAction.h ./source/actions/ImprovedGaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ImprovedGaugeAction.o ./source/actions/ImprovedGaugeAction.cpp

./build/GaugeForce.o: ./source/hmc_forces/GaugeForce.h ./source/hmc_forces/GaugeForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GaugeForce.o ./source/hmc_forces/GaugeForce.cpp

./build/GaugeAction.o: ./source/actions/GaugeAction.h ./source/actions/GaugeAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GaugeAction.o ./source/actions/GaugeAction.cpp

./build/Integrate.o: ./source/hmc_integrators/Integrate.h ./source/hmc_integrators/Integrate.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Integrate.o ./source/hmc_integrators/Integrate.cpp

./build/LeapFrog.o: ./source/hmc_integrators/LeapFrog.h ./source/hmc_integrators/LeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LeapFrog.o ./source/hmc_integrators/LeapFrog.cpp

./build/FourthOrderLeapFrog.o: ./source/hmc_integrators/FourthOrderLeapFrog.h ./source/hmc_integrators/FourthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/FourthOrderLeapFrog.o ./source/hmc_integrators/FourthOrderLeapFrog.cpp

./build/SixthOrderLeapFrog.o: ./source/hmc_integrators/SixthOrderLeapFrog.h ./source/hmc_integrators/SixthOrderLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SixthOrderLeapFrog.o ./source/hmc_integrators/SixthOrderLeapFrog.cpp

./build/OmelyanLeapFrog.o: ./source/hmc_integrators/OmelyanLeapFrog.h ./source/hmc_integrators/OmelyanLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/OmelyanLeapFrog.o ./source/hmc_integrators/OmelyanLeapFrog.cpp

./build/FourthOmelyanLeapFrog.o: ./source/hmc_integrators/FourthOmelyanLeapFrog.h ./source/hmc_integrators/FourthOmelyanLeapFrog.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/FourthOmelyanLeapFrog.o ./source/hmc_integrators/FourthOmelyanLeapFrog.cpp

./build/Energy.o: ./source/actions/Energy.h ./source/actions/Energy.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Energy.o ./source/actions/Energy.cpp

./build/Force.o: ./source/hmc_forces/Force.h ./source/hmc_forces/Force.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Force.o ./source/hmc_forces/Force.cpp

./build/HMCUpdater.o: ./source/hmc_updaters/HMCUpdater.h ./source/hmc_updaters/HMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/HMCUpdater.o ./source/hmc_updaters/HMCUpdater.cpp

./build/FermionHMCUpdater.o: ./source/hmc_updaters/FermionHMCUpdater.h ./source/hmc_updaters/FermionHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/FermionHMCUpdater.o ./source/hmc_updaters/FermionHMCUpdater.cpp

./build/FermionicAction.o: ./source/actions/FermionicAction.h ./source/actions/FermionicAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/FermionicAction.o ./source/actions/FermionicAction.cpp

./build/Polynomial.o: ./source/dirac_functions/Polynomial.h ./source/dirac_functions/Polynomial.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Polynomial.o ./source/dirac_functions/Polynomial.cpp

./build/ChebyshevRecursion.o: ./source/dirac_functions/ChebyshevRecursion.h ./source/dirac_functions/ChebyshevRecursion.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ChebyshevRecursion.o ./source/dirac_functions/ChebyshevRecursion.cpp

./build/RationalApproximation.o: ./source/dirac_functions/RationalApproximation.h ./source/dirac_functions/RationalApproximation.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/RationalApproximation.o ./source/dirac_functions/RationalApproximation.cpp

./build/GaugeFixing.o: ./source/gauge_fixing/GaugeFixing.h ./source/gauge_fixing/GaugeFixing.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GaugeFixing.o ./source/gauge_fixing/GaugeFixing.cpp

./build/LandauGaugeFixing.o: ./source/gauge_fixing/LandauGaugeFixing.h ./source/gauge_fixing/LandauGaugeFixing.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LandauGaugeFixing.o ./source/gauge_fixing/LandauGaugeFixing.cpp

./build/LandauGluonPropagator.o: ./source/gauge_fixing/LandauGluonPropagator.h ./source/gauge_fixing/LandauGluonPropagator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LandauGluonPropagator.o ./source/gauge_fixing/LandauGluonPropagator.cpp

./build/LandauGhostPropagator.o: ./source/gauge_fixing/LandauGhostPropagator.h ./source/gauge_fixing/LandauGhostPropagator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/LandauGhostPropagator.o ./source/gauge_fixing/LandauGhostPropagator.cpp

./build/MaximalAbelianGaugeFixing.o: ./source/gauge_fixing/MaximalAbelianGaugeFixing.h ./source/gauge_fixing/MaximalAbelianGaugeFixing.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MaximalAbelianGaugeFixing.o ./source/gauge_fixing/MaximalAbelianGaugeFixing.cpp

./build/MaximalAbelianProjection.o: ./source/gauge_fixing/MaximalAbelianProjection.h ./source/gauge_fixing/MaximalAbelianProjection.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MaximalAbelianProjection.o ./source/gauge_fixing/MaximalAbelianProjection.cpp

./build/ReUnit.o: ./source/utils/ReUnit.h ./source/utils/ReUnit.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ReUnit.o ./source/utils/ReUnit.cpp

./build/RandomGaugeTransformation.o: ./source/utils/RandomGaugeTransformation.h ./source/utils/RandomGaugeTransformation.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/RandomGaugeTransformation.o ./source/utils/RandomGaugeTransformation.cpp

./build/StoutSmearing.o: ./source/utils/StoutSmearing.h ./source/utils/StoutSmearing.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/StoutSmearing.o ./source/utils/StoutSmearing.cpp

./build/Gamma.o: ./source/utils/Gamma.h ./source/utils/Gamma.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Gamma.o ./source/utils/Gamma.cpp

./build/PureGaugeUpdater.o: ./source/pure_gauge/PureGaugeUpdater.h ./source/pure_gauge/PureGaugeUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PureGaugeUpdater.o ./source/pure_gauge/PureGaugeUpdater.cpp

./build/PureGaugeOverrelaxation.o: ./source/pure_gauge/PureGaugeOverrelaxation.h ./source/pure_gauge/PureGaugeOverrelaxation.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PureGaugeOverrelaxation.o ./source/pure_gauge/PureGaugeOverrelaxation.cpp

./build/PureGaugeHMCUpdater.o: ./source/hmc_updaters/PureGaugeHMCUpdater.h ./source/hmc_updaters/PureGaugeHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PureGaugeHMCUpdater.o ./source/hmc_updaters/PureGaugeHMCUpdater.cpp

./build/Checkerboard.o: ./source/pure_gauge/Checkerboard.h ./source/pure_gauge/Checkerboard.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Checkerboard.o ./source/pure_gauge/Checkerboard.cpp

./build/RandomSeed.o: ./source/utils/RandomSeed.h ./source/utils/RandomSeed.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/RandomSeed.o ./source/utils/RandomSeed.cpp

./build/Plaquette.o: ./source/wilson_loops/Plaquette.h ./source/wilson_loops/Plaquette.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Plaquette.o ./source/wilson_loops/Plaquette.cpp

./build/PolyakovLoop.o: ./source/polyakov_loops/PolyakovLoop.h ./source/polyakov_loops/PolyakovLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PolyakovLoop.o ./source/polyakov_loops/PolyakovLoop.cpp

./build/PolyakovLoopEigenvalues.o: ./source/polyakov_loops/PolyakovLoopEigenvalues.h ./source/polyakov_loops/PolyakovLoopEigenvalues.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PolyakovLoopEigenvalues.o ./source/polyakov_loops/PolyakovLoopEigenvalues.cpp

./build/PolyakovLoopCorrelator.o: ./source/polyakov_loops/PolyakovLoopCorrelator.h ./source/polyakov_loops/PolyakovLoopCorrelator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PolyakovLoopCorrelator.o ./source/polyakov_loops/PolyakovLoopCorrelator.cpp

./build/AdjointPolyakovLoop.o: ./source/polyakov_loops/AdjointPolyakovLoop.h ./source/polyakov_loops/AdjointPolyakovLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/AdjointPolyakovLoop.o ./source/polyakov_loops/AdjointPolyakovLoop.cpp

./build/WilsonLoop.o: ./source/wilson_loops/WilsonLoop.h ./source/wilson_loops/WilsonLoop.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/WilsonLoop.o ./source/wilson_loops/WilsonLoop.cpp

./build/PureGaugeWilsonLoops.o: ./source/pure_gauge/PureGaugeWilsonLoops.h ./source/pure_gauge/PureGaugeWilsonLoops.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/PureGaugeWilsonLoops.o ./source/pure_gauge/PureGaugeWilsonLoops.cpp 

./build/StartGaugeConfiguration.o: ./source/starters/StartGaugeConfiguration.h ./source/starters/StartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/StartGaugeConfiguration.o ./source/starters/StartGaugeConfiguration.cpp

./build/ReadStartGaugeConfiguration.o: ./source/starters/ReadStartGaugeConfiguration.h ./source/starters/ReadStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ReadStartGaugeConfiguration.o ./source/starters/ReadStartGaugeConfiguration.cpp

./build/HotStartGaugeConfiguration.o: ./source/starters/HotStartGaugeConfiguration.h ./source/starters/HotStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/HotStartGaugeConfiguration.o ./source/starters/HotStartGaugeConfiguration.cpp

./build/ColdStartGaugeConfiguration.o: ./source/starters/ColdStartGaugeConfiguration.h ./source/starters/ColdStartGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ColdStartGaugeConfiguration.o ./source/starters/ColdStartGaugeConfiguration.cpp

./build/StochasticEstimator.o: ./source/fermion_measurements/StochasticEstimator.h ./source/fermion_measurements/StochasticEstimator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/StochasticEstimator.o ./source/fermion_measurements/StochasticEstimator.cpp

./build/MesonCorrelator.o: ./source/correlators/MesonCorrelator.h ./source/correlators/MesonCorrelator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MesonCorrelator.o ./source/correlators/MesonCorrelator.cpp

./build/OverlapChiralRotation.o: ./source/fermion_measurements/OverlapChiralRotation.h ./source/fermion_measurements/OverlapChiralRotation.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/OverlapChiralRotation.o ./source/fermion_measurements/OverlapChiralRotation.cpp

./build/ChiralCondensate.o: ./source/fermion_measurements/ChiralCondensate.h ./source/fermion_measurements/ChiralCondensate.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ChiralCondensate.o ./source/fermion_measurements/ChiralCondensate.cpp

./build/NPRVertex.o: ./source/fermion_measurements/NPRVertex.h ./source/fermion_measurements/NPRVertex.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/NPRVertex.o ./source/fermion_measurements/NPRVertex.cpp
	
./build/GluinoGlue.o: ./source/correlators/GluinoGlue.h ./source/correlators/GluinoGlue.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GluinoGlue.o ./source/correlators/GluinoGlue.cpp

./build/FermionForce.o: ./source/hmc_forces/FermionForce.h ./source/hmc_forces/FermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/FermionForce.o ./source/hmc_forces/FermionForce.cpp

./build/SmearingForce.o: ./source/hmc_forces/SmearingForce.h ./source/hmc_forces/SmearingForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SmearingForce.o ./source/hmc_forces/SmearingForce.cpp

./build/DiracWilsonFermionForce.o: ./source/hmc_forces/DiracWilsonFermionForce.h ./source/hmc_forces/DiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/DiracWilsonFermionForce.o ./source/hmc_forces/DiracWilsonFermionForce.cpp

./build/OverlapFermionForce.o: ./source/hmc_forces/OverlapFermionForce.h ./source/hmc_forces/OverlapFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/OverlapFermionForce.o ./source/hmc_forces/OverlapFermionForce.cpp

./build/BlockDiracWilsonFermionForce.o: ./source/hmc_forces/BlockDiracWilsonFermionForce.h ./source/hmc_forces/BlockDiracWilsonFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/BlockDiracWilsonFermionForce.o ./source/hmc_forces/BlockDiracWilsonFermionForce.cpp

./build/ImprovedFermionForce.o: ./source/hmc_forces/ImprovedFermionForce.h ./source/hmc_forces/ImprovedFermionForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ImprovedFermionForce.o ./source/hmc_forces/ImprovedFermionForce.cpp

./build/TestForce.o: ./source/hmc_forces/TestForce.h ./source/hmc_forces/TestForce.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TestForce.o ./source/hmc_forces/TestForce.cpp

./build/ScalarAction.o: ./source/actions/ScalarAction.h ./source/actions/ScalarAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ScalarAction.o ./source/actions/ScalarAction.cpp

./build/MultiScalarAction.o: ./source/actions/MultiScalarAction.h ./source/actions/MultiScalarAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiScalarAction.o ./source/actions/MultiScalarAction.cpp

./build/TwoFlavorFermionAction.o: ./source/actions/TwoFlavorFermionAction.h ./source/actions/TwoFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TwoFlavorFermionAction.o ./source/actions/TwoFlavorFermionAction.cpp

./build/TwoFlavorQCDAction.o: ./source/actions/TwoFlavorQCDAction.h ./source/actions/TwoFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TwoFlavorQCDAction.o ./source/actions/TwoFlavorQCDAction.cpp

./build/TwoFlavorHMCUpdater.o: ./source/hmc_updaters/TwoFlavorHMCUpdater.h ./source/hmc_updaters/TwoFlavorHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TwoFlavorHMCUpdater.o ./source/hmc_updaters/TwoFlavorHMCUpdater.cpp

./build/ScalarFermionHMCUpdater.o: ./source/hmc_updaters/ScalarFermionHMCUpdater.h ./source/hmc_updaters/ScalarFermionHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ScalarFermionHMCUpdater.o ./source/hmc_updaters/ScalarFermionHMCUpdater.cpp

./build/RandomScalarUpdater.o: ./source/scalar_updaters/RandomScalarUpdater.h ./source/scalar_updaters/RandomScalarUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/RandomScalarUpdater.o ./source/scalar_updaters/RandomScalarUpdater.cpp

./build/MeanScalarField.o: ./source/scalar_measurements/MeanScalarField.cpp ./source/scalar_measurements/MeanScalarField.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MeanScalarField.o ./source/scalar_measurements/MeanScalarField.cpp

./build/MetropolisScalarUpdater.o: ./source/scalar_updaters/MetropolisScalarUpdater.h ./source/scalar_updaters/MetropolisScalarUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MetropolisScalarUpdater.o ./source/scalar_updaters/MetropolisScalarUpdater.cpp

./build/NFlavorFermionAction.o: ./source/actions/NFlavorFermionAction.h ./source/actions/NFlavorFermionAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/NFlavorFermionAction.o ./source/actions/NFlavorFermionAction.cpp

./build/NFlavorQCDAction.o: ./source/actions/NFlavorQCDAction.h ./source/actions/NFlavorQCDAction.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/NFlavorQCDAction.o ./source/actions/NFlavorQCDAction.cpp

./build/HiggsGaugeHMCUpdater.o : ./source/hmc_updaters/HiggsGaugeHMCUpdater.h ./source/hmc_updaters/HiggsGaugeHMCUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/HiggsGaugeHMCUpdater.o ./source/hmc_updaters/HiggsGaugeHMCUpdater.cpp

./build/MultiStepNFlavorUpdater.o : ./source/hmc_updaters/MultiStepNFlavorUpdater.h ./source/hmc_updaters/MultiStepNFlavorUpdater.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiStepNFlavorUpdater.o ./source/hmc_updaters/MultiStepNFlavorUpdater.cpp

./build/DiracEigenSolver.o: ./source/fermion_measurements/DiracEigenSolver.h ./source/fermion_measurements/DiracEigenSolver.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/DiracEigenSolver.o ./source/fermion_measurements/DiracEigenSolver.cpp

./build/SingletOperators.o: ./source/fermion_measurements/SingletOperators.h ./source/fermion_measurements/SingletOperators.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SingletOperators.o ./source/fermion_measurements/SingletOperators.cpp

./build/XSpaceCorrelators.o: ./source/fermion_measurements/XSpaceCorrelators.h ./source/fermion_measurements/XSpaceCorrelators.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/XSpaceCorrelators.o ./source/fermion_measurements/XSpaceCorrelators.cpp

./build/Eigenvalues.o: ./source/fermion_measurements/Eigenvalues.h ./source/fermion_measurements/Eigenvalues.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Eigenvalues.o ./source/fermion_measurements/Eigenvalues.cpp

./build/GlobalOutput.o: ./source/io/GlobalOutput.h ./source/io/GlobalOutput.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GlobalOutput.o ./source/io/GlobalOutput.cpp

./build/OutputSweep.o: ./source/io/OutputSweep.h ./source/io/OutputSweep.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/OutputSweep.o ./source/io/OutputSweep.cpp

./build/Glueball.o: ./source/correlators/Glueball.h ./source/correlators/Glueball.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/Glueball.o ./source/correlators/Glueball.cpp

./build/TestCommunication.o: ./source/tests/TestCommunication.h ./source/tests/TestCommunication.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TestCommunication.o ./source/tests/TestCommunication.cpp

./build/TestLinearAlgebra.o: ./source/tests/TestLinearAlgebra.h ./source/tests/TestLinearAlgebra.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TestLinearAlgebra.o ./source/tests/TestLinearAlgebra.cpp

./build/TestSpeedDiracOperators.o: ./source/tests/TestSpeedDiracOperators.h ./source/tests/TestSpeedDiracOperators.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/TestSpeedDiracOperators.o ./source/tests/TestSpeedDiracOperators.cpp

./build/StorageParameters.o: ./source/program_options/StorageParameters.h ./source/program_options/StorageParameters.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/StorageParameters.o ./source/program_options/StorageParameters.cpp

./build/ReadGaugeConfiguration.o: ./source/io/ReadGaugeConfiguration.h ./source/io/ReadGaugeConfiguration.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/ReadGaugeConfiguration.o ./source/io/ReadGaugeConfiguration.cpp
	
./build/WilsonFlow.o: ./source/wilson_flow/WilsonFlow.h ./source/wilson_flow/WilsonFlow.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/WilsonFlow.o ./source/wilson_flow/WilsonFlow.cpp

./build/GaugeEnergy.o: ./source/actions/GaugeEnergy.h ./source/actions/GaugeEnergy.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/GaugeEnergy.o ./source/actions/GaugeEnergy.cpp

./build/MultiGridDiracOperator.o: ./source/multigrid_solver/MultiGridDiracOperator.h ./source/multigrid_solver/MultiGridDiracOperator.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridDiracOperator.o ./source/multigrid_solver/MultiGridDiracOperator.cpp

./build/SAPPreconditioner.o: ./source/dirac_operators/SAPPreconditioner.h ./source/dirac_operators/SAPPreconditioner.cpp
	$(CPP) $(CPPFLAGS) -c -o ./build/SAPPreconditioner.o ./source/dirac_operators/SAPPreconditioner.cpp
	
./build/BlockBasis.o: ./source/multigrid/BlockBasis.cpp ./source/multigrid/BlockBasis.h
	$(CPP) $(CPPFLAGS) -c -o ./build/BlockBasis.o ./source/multigrid/BlockBasis.cpp

./build/MultiGridBiConjugateGradient.o: ./source/multigrid/MultiGridBiConjugateGradient.cpp ./source/multigrid/MultiGridBiConjugateGradient.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridBiConjugateGradient.o ./source/multigrid/MultiGridBiConjugateGradient.cpp

./build/MultiGridConjugateGradient.o: ./source/multigrid/MultiGridConjugateGradient.cpp ./source/multigrid/MultiGridConjugateGradient.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridConjugateGradient.o ./source/multigrid/MultiGridConjugateGradient.cpp

./build/MultiGridOperator.o: ./source/multigrid/MultiGridOperator.cpp ./source/multigrid/MultiGridOperator.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridOperator.o ./source/multigrid/MultiGridOperator.cpp

./build/MultiGridProjector.o: ./source/multigrid/MultiGridProjector.cpp ./source/multigrid/MultiGridProjector.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridProjector.o ./source/multigrid/MultiGridProjector.cpp

./build/MultiGridSolver.o: ./source/multigrid/MultiGridSolver.cpp ./source/multigrid/MultiGridSolver.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridSolver.o ./source/multigrid/MultiGridSolver.cpp

./build/MultiGridVectorLayout.o: ./source/multigrid/MultiGridVectorLayout.cpp ./source/multigrid/MultiGridVectorLayout.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridVectorLayout.o ./source/multigrid/MultiGridVectorLayout.cpp

./build/MultiGridStochasticEstimator.o: ./source/multigrid/MultiGridStochasticEstimator.cpp ./source/multigrid/MultiGridStochasticEstimator.h
	$(CPP) $(CPPFLAGS) -c -o ./build/MultiGridStochasticEstimator.o ./source/multigrid/MultiGridStochasticEstimator.cpp


