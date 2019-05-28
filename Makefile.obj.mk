OBJECTS  = ReducedStencil.o StandardStencil.o ExtendedStencil.o LocalLayout.o \
			AlgebraUtils.o \
			BiConjugateGradient.o DeflationInverter.o ConjugateGradient.o MultishiftSolver.o ChronologicalMultishiftSolver.o MMMRMultishiftSolver.o MEMultishiftSolver.o MultiGridMEMultishiftSolver.o GMRESR.o PreconditionedBiCGStab.o \
			AdjointScalarAction.o FundamentalScalarAction.o ScalarAction.o MultiScalarAction.o \
			DiracOperator.o Propagator.o BasicDiracWilsonOperator.o BasicSquareDiracWilsonOperator.o DiracWilsonOperator.o SquareDiracWilsonOperator.o BlockDiracWilsonOperator.o BlockImprovedDiracWilsonOperator.o BlockDiracOperator.o ComplementBlockDiracOperator.o OverlapOperator.o SquareOverlapOperator.o ExactOverlapOperator.o SquareComplementBlockDiracWilsonOperator.o SquareComplementBlockDiracOperator.o SquareBlockDiracWilsonOperator.o ImprovedDiracWilsonOperator.o SquareImprovedDiracWilsonOperator.o SquareTwistedDiracOperator.o TwistedDiracOperator.o SAPPreconditioner.o HoppingOperator.o GammaOperators.o EvenOddImprovedDiracWilsonOperator.o SquareEvenOddImprovedDiracWilsonOperator.o \
			BlockBasis.o MultiGridBiConjugateGradient.o MultiGridConjugateGradient.o MultiGridOperator.o MultiGridProjector.o MultiGridSolver.o MultiGridVectorLayout.o MultiGridStochasticEstimator.o \
			Polynomial.o RationalApproximation.o ChebyshevRecursion.o \
			Integrate.o LeapFrog.o FourthOrderLeapFrog.o SixthOrderLeapFrog.o OmelyanLeapFrog.o FourthOmelyanLeapFrog.o Energy.o Force.o \
			HMCUpdater.o FermionHMCUpdater.o \
			ScalarFermionHMCUpdater.o RandomScalarUpdater.o AdjointMetropolisScalarUpdater.o MeanScalarField.o FundamentalMetropolisScalarUpdater.o HiggsGaugeHMCUpdater.o \
			FermionicAction.o \
			GaugeForce.o GaugeAction.o WilsonGaugeAction.o ImprovedGaugeAction.o \
			ReUnit.o StoutSmearing.o Gamma.o \
			RandomSeed.o \
			GaugeFixing.o LandauGaugeFixing.o MaximalAbelianGaugeFixing.o MaximalAbelianProjection.o LandauGluonPropagator.o LandauGhostPropagator.o \
			Glueball.o \
			Plaquette.o PolyakovLoop.o PolyakovLoopEigenvalues.o PolyakovLoopCorrelator.o AdjointPolyakovLoop.o WilsonLoop.o GaugeEnergy.o \
			GlobalOutput.o OutputSweep.o \
			FermionForce.o DiracWilsonFermionForce.o BlockDiracWilsonFermionForce.o ImprovedFermionForce.o TestForce.o SmearingForce.o OverlapFermionForce.o \
			StochasticEstimator.o MesonCorrelator.o ChiralCondensate.o SingletOperators.o GluinoGlue.o NPRVertex.o XSpaceCorrelators.o OverlapResidualMass.o \
			PureGaugeUpdater.o PureGaugeOverrelaxation.o PureGaugeHMCUpdater.o Checkerboard.o PureGaugeWilsonLoops.o \
			TwoFlavorFermionAction.o TwoFlavorQCDAction.o TwoFlavorHMCUpdater.o \
			NFlavorFermionAction.o NFlavorQCDAction.o MultiStepNFlavorQCDUpdater.o \
			DiracEigenSolver.o Eigenvalues.o \
			TestCommunication.o TestLinearAlgebra.o TestSpeedDiracOperators.o \
			StartGaugeConfiguration.o ReadStartGaugeConfiguration.o HotStartGaugeConfiguration.o ColdStartGaugeConfiguration.o \
			ReadGaugeConfiguration.o \
			LatticeSweep.o Simulation.o \
			StorageParameters.o \
			WilsonFlow.o \
			Environment.o updater.o \

