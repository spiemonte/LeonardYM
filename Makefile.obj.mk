OBJECTS  = ReducedStencil.o StandardStencil.o ExtendedStencil.o LocalLayout.o \
			AlgebraUtils.o \
			BiConjugateGradient.o DeflationInverter.o ConjugateGradient.o MultishiftSolver.o ChronologicalMultishiftSolver.o MMMRMultishiftSolver.o MEMultishiftSolver.o \
			DiracOperator.o BasicDiracWilsonOperator.o BasicSquareDiracWilsonOperator.o DiracWilsonOperator.o SquareDiracWilsonOperator.o BlockDiracWilsonOperator.o BlockDiracOperator.o ComplementBlockDiracWilsonOperator.o SquareComplementBlockDiracWilsonOperator.o SquareBlockDiracWilsonOperator.o ImprovedDiracWilsonOperator.o SquareImprovedDiracWilsonOperator.o SquareTwistedDiracOperator.o \
			Polynomial.o RationalApproximation.o \
			Integrate.o LeapFrog.o FourthOrderLeapFrog.o SixthOrderLeapFrog.o OmelyanLeapFrog.o FourthOmelyanLeapFrog.o Energy.o Force.o \
			HMCUpdater.o FermionHMCUpdater.o \
			BandAction.o BandTwoFlavorUpdater.o \
			FermionicAction.o \
			GaugeForce.o GaugeAction.o WilsonGaugeAction.o ImprovedGaugeAction.o \
			ReUnit.o StoutSmearing.o Gamma.o \
			RandomSeed.o \
			Glueball.o \
			Plaquette.o PolyakovLoop.o AdjointPolyakovLoop.o  WilsonLoop.o GaugeEnergy.o \
			GlobalOutput.o OutputSweep.o \
			FermionForce.o DiracWilsonFermionForce.o BlockDiracWilsonFermionForce.o ImprovedFermionForce.o TestForce.o \
			StochasticEstimator.o MesonCorrelator.o ChiralCondensate.o GluinoGlue.o \
			PureGaugeUpdater.o PureGaugeOverrelaxation.o PureGaugeHMCUpdater.o Checkerboard.o \
			TwoFlavorFermionAction.o TwoFlavorQCDAction.o TwoFlavorHMCUpdater.o \
			NFlavorFermionAction.o NFlavorQCDAction.o NFlavorQCDUpdater.o MultiStepNFlavorQCDUpdater.o TwistedMultiStepNFlavorQCDUpdater.o NFlavorBlockUpdater.o \
			DiracEigenSolver.o Eigenvalues.o \
			TestCommunication.o TestLinearAlgebra.o \
			StartGaugeConfiguration.o ReadStartGaugeConfiguration.o HotStartGaugeConfiguration.o ColdStartGaugeConfiguration.o \
			ReadGaugeConfiguration.o \
			LatticeSweep.o Simulation.o \
			StorageParameters.o \
			WilsonFlow.o \
			updater.o \
			#MultiGridDiracOperator.o MultiGridSolver.o \
