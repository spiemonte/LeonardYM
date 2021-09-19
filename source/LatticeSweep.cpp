#include "LatticeSweep.h"
#include "pure_gauge/PureGaugeUpdater.h"
#include "wilson_loops/Plaquette.h"
#include "utils/FromString.h"
#include "hmc_updaters/PureGaugeHMCUpdater.h"
#include "hmc_updaters/TwoFlavorHMCUpdater.h"
#include "hmc_updaters/ScalarFermionHMCUpdater.h"
#include "pure_gauge/PureGaugeOverrelaxation.h"
#include "io/OutputSweep.h"
#include "tests/TestLinearAlgebra.h"
#include "tests/TestSpeedDiracOperators.h"
#include "fermion_measurements/Eigenvalues.h"
#include "correlators/MesonCorrelator.h"
#include "fermion_measurements/ChiralCondensate.h"
#include "fermion_measurements/OverlapChiralRotation.h"
#include "fermion_measurements/SingletOperators.h"
#include "fermion_measurements/XSpaceCorrelators.h"
#include "fermion_measurements/NPRVertex.h"
#include "scalar_measurements/MeanScalarField.h"
#include "polyakov_loops/PolyakovLoop.h"
#include "correlators/Glueball.h"
#include "utils/ReUnit.h"
#include "gauge_fixing/LandauGaugeFixing.h"
#include "gauge_fixing/MaximalAbelianGaugeFixing.h"
#include "gauge_fixing/MaximalAbelianProjection.h"
#include "gauge_fixing/LandauGhostPropagator.h"
#include "gauge_fixing/LandauGluonPropagator.h"
#include "pure_gauge/PureGaugeWilsonLoops.h"
#include "tests/TestCommunication.h"
#include "correlators/GluinoGlue.h"
#include "wilson_loops/WilsonLoop.h"
#include "io/ReadGaugeConfiguration.h"
#include "wilson_flow/WilsonFlow.h"
#include "actions/GaugeEnergy.h"
#include "polyakov_loops/AdjointPolyakovLoop.h"
#include "polyakov_loops/PolyakovLoopCorrelator.h"
#include "polyakov_loops/PolyakovLoopEigenvalues.h"
#include "hmc_updaters/MultiStepNFlavorUpdater.h"
#include "scalar_updaters/RandomScalarUpdater.h"
#include "scalar_updaters/MetropolisScalarUpdater.h"
#include "hmc_updaters/HiggsGaugeHMCUpdater.h"
#include "utils/RandomGaugeTransformation.h"

namespace Update {

	LatticeSweep::LatticeSweep() { }

	LatticeSweep::LatticeSweep(unsigned int _numberTimes, unsigned int _sweepToJump) : numberTimes(_numberTimes), sweepToJump(_sweepToJump) { }

	LatticeSweep::~LatticeSweep() { }

	LatticeSweep* LatticeSweep::getInstance(const std::string& name) {
	//Factory for the Lattice sweep
		if (name == "PureGaugeCM") {
			return new PureGaugeUpdater();
		} else if (name == "Plaquette") {
			return new Plaquette();
		} else if (name == "PureGaugeHMC") {
			return new PureGaugeHMCUpdater();
		} else if (name == "TwoFlavorQCD") {
			return new TwoFlavorHMCUpdater();
		} else if (name == "PureGaugeOverrelaxation") {
			return new PureGaugeOverrelaxation();
		} else if (name == "Output") {
			return new OutputSweep();
		} else if (name == "TestLinearAlgebra") {
			return new TestLinearAlgebra();
		} else if (name == "TestSpeedDiracOperators") {
			return new TestSpeedDiracOperators();
		} else if (name == "Eigenvalues") {
			return new Eigenvalues();
		} else if (name == "MesonCorrelator") {
			return new MesonCorrelator();
		} else if (name == "OverlapChiralRotation") {
			return new OverlapChiralRotation();
		} else if (name == "ChiralCondensate") {
			return new ChiralCondensate();
		} else if (name == "PolyakovLoop") {
			return new PolyakovLoop();
		} else if (name == "PolyakovLoopEigenvalues") {
			return new PolyakovLoopEigenvalues();
		} else if (name == "PolyakovLoopCorrelator") {
			return new PolyakovLoopCorrelator();
		} else if (name == "AdjointPolyakovLoop") {
			return new AdjointPolyakovLoop();
		} else if (name == "Glueball") {
			return new Glueball();
		} else if (name == "GluinoGlue") {
			return new GluinoGlue();
		} else if (name == "ReUnit") {
			return new ReUnit();
		} else if (name == "PureGaugeWilsonLoops") {
			return new PureGaugeWilsonLoops();
		} else if (name == "WilsonLoop") {
			return new WilsonLoop();
		} else if (name == "TestCommunication") {
			return new TestCommunication();
		} else if (name == "ReadGaugeConfiguration") {
			return new ReadGaugeConfiguration();
		} else if (name == "WilsonFlow") {
			return new WilsonFlow();
		} else if (name == "GaugeEnergy") {
			return new GaugeEnergy();
		} else if (name == "MultiStepNFlavor") {
			return new MultiStepNFlavorUpdater();
		} else if (name == "SingletOperators") {
			return new SingletOperators();
		} else if (name == "XSpaceCorrelators") {
			return new XSpaceCorrelators();
		} else if (name == "NPRVertex") {
			return new NPRVertex();
		} else if (name == "LandauGaugeFixing") {
			return new LandauGaugeFixing();
		} else if (name == "MaximalAbelianGaugeFixing") {
			return new MaximalAbelianGaugeFixing();
		} else if (name == "MaximalAbelianProjection") {
			return new MaximalAbelianProjection();
		} else if (name == "LandauGluonPropagator") {
			return new LandauGluonPropagator();
		} else if (name == "LandauGhostPropagator") {
			return new LandauGhostPropagator();
		} else if (name == "ScalarFermionHMC") {
			return new ScalarFermionHMCUpdater();
		} else if (name == "HiggsGaugeHMC") {
			return new HiggsGaugeHMCUpdater();
		} else if (name == "RandomScalarInitializer") {
			return new RandomScalarUpdater();
		} else if (name == "MCScalar") {
			return new MetropolisScalarUpdater();
		} else if (name == "MeanScalarField") {
			return new MeanScalarField();
		} else if (name == "RandomGaugeTransformation") {
			return new RandomGaugeTransformation();
		}
		else {
			if (isOutputProcess()) std::cout << "Unknown name sweep: " << name << std::endl;
			exit(1);
		}
	}

	void LatticeSweep::call(environment_t& environment) {
	//Execute the sweeps in the order given
		if ((environment.sweep % sweepToJump) == 0) {
			environment.iteration = 0;
			for (unsigned int i = 0; i < numberTimes; ++i) {
				this->execute(environment);
				++environment.iteration;
			}
		}
	}

	LatticeSweep* LatticeSweep::read(const std::string& toRead) {
		std::string::const_iterator index = toRead.begin();
		while (index != toRead.end() && *index != '{') {
			++index;
		}
		if (index != toRead.end()) ++index;
		else {
			std::cout << "Bad lattice sweep vector format, did you forget to open {?" << std::endl;
			exit(17);
		}

		std::string nameSweep;
		while (index != toRead.end() && *index != ',') {
			nameSweep.push_back(*index);
			++index;
		}
		if (index != toRead.end()) ++index;
		else {
			std::cout << "Bad lattice sweep vector format, did you forget to put some comma?" << std::endl;
			exit(17);
		}

		std::string steps;
		while (index != toRead.end() && *index != ',') {
			steps.push_back(*index);
			++index;
		}
		if (index != toRead.end()) ++index;
		else {
			std::cout << "Bad lattice sweep vector format, did you forget to put some comma?" << std::endl;
			exit(17);
		}

		std::string numberCalls;
		while (index != toRead.end() && *index != '}') {
			numberCalls.push_back(*index);
			++index;
		}
		if (index != toRead.end()) ++index;
		else {
			std::cout << "Bad lattice sweep vector format, did you forget to close }?" << std::endl;
			exit(17);
		}

		LatticeSweep* sweep = LatticeSweep::getInstance(nameSweep);
		sweep->setSweepToJump(fromString<int>(steps));
		sweep->setNumberTimes(fromString<int>(numberCalls));
		return sweep;
	}

	void LatticeSweep::setNumberTimes(unsigned int _numberTimes) {
		numberTimes = _numberTimes;
	}

	void LatticeSweep::setSweepToJump(unsigned int _sweepToJump) {
		sweepToJump = _sweepToJump;
	}

	void LatticeSweep::addParameters(std::map<std::string, Option>& desc) {
		PureGaugeUpdater::registerParameters(desc);
		Plaquette::registerParameters(desc);
		PureGaugeHMCUpdater::registerParameters(desc);
		TwoFlavorHMCUpdater::registerParameters(desc);
		PureGaugeOverrelaxation::registerParameters(desc);
		OutputSweep::registerParameters(desc);
		TestLinearAlgebra::registerParameters(desc);
		TestSpeedDiracOperators::registerParameters(desc);
		Eigenvalues::registerParameters(desc);
		MesonCorrelator::registerParameters(desc);
		ChiralCondensate::registerParameters(desc);
		OverlapChiralRotation::registerParameters(desc);
		PolyakovLoopEigenvalues::registerParameters(desc);
		PolyakovLoopCorrelator::registerParameters(desc);
		AdjointPolyakovLoop::registerParameters(desc);
		Glueball::registerParameters(desc);
		GluinoGlue::registerParameters(desc);
		ReUnit::registerParameters(desc);
		PureGaugeWilsonLoops::registerParameters(desc);
		WilsonLoop::registerParameters(desc);
		TestCommunication::registerParameters(desc);
		ReadGaugeConfiguration::registerParameters(desc);
		WilsonFlow::registerParameters(desc);
		GaugeEnergy::registerParameters(desc);
		MultiStepNFlavorUpdater::registerParameters(desc);
		SingletOperators::registerParameters(desc);
		XSpaceCorrelators::registerParameters(desc);
		NPRVertex::registerParameters(desc);
		LandauGaugeFixing::registerParameters(desc);
		MaximalAbelianGaugeFixing::registerParameters(desc);
		MaximalAbelianProjection::registerParameters(desc);
		LandauGluonPropagator::registerParameters(desc);
		LandauGhostPropagator::registerParameters(desc);
		ScalarFermionHMCUpdater::registerParameters(desc);
		HiggsGaugeHMCUpdater::registerParameters(desc);
		MetropolisScalarUpdater::registerParameters(desc);
		RandomGaugeTransformation::registerParameters(desc);
	}

	void LatticeSweep::registerParameters(std::map<std::string, Option>&) {
	}

	void LatticeSweep::printSweepsName() {
		if (isOutputProcess()) {
			std::cout << "List of the possible Sweeps name:" << std::endl;
			std::cout
			<< "PureGaugeCM" << std::endl
			<<  "Plaquette" << std::endl
			<<  "PureGaugeHMC" << std::endl
			<<  "TwoFlavorQCD" << std::endl
			<<  "PureGaugeOverrelaxation" << std::endl
			<<  "Output" << std::endl
			<<  "TestLinearAlgebra" << std::endl
			<<  "TestSpeedDiracOperators" << std::endl
			<<  "Eigenvalues" << std::endl
			<<  "MesonCorrelator" << std::endl
			<<  "OverlapChiralRotation" << std::endl
			<<  "ChiralCondensate" << std::endl
			<<  "PolyakovLoop" << std::endl
			<<  "PolyakovLoopEigenvalues" << std::endl
			<<  "PolyakovLoopCorrelator" << std::endl
			<<  "AdjointPolyakovLoop" << std::endl
			<<  "Glueball" << std::endl
			<<  "GluinoGlue" << std::endl
			<<  "ReUnit" << std::endl
			<<  "PureGaugeWilsonLoops" << std::endl
			<<  "WilsonLoop" << std::endl
			<<  "TestCommunication" << std::endl
			<<  "ReadGaugeConfiguration" << std::endl
			<<  "WilsonFlow" << std::endl
			<<  "GaugeEnergy" << std::endl
			<<  "MultiStepNFlavor" << std::endl
			<<  "SingletOperators" << std::endl
			<<  "XSpaceCorrelators" << std::endl
			<<  "NPRVertex" << std::endl
			<<  "LandauGaugeFixing" << std::endl
			<<  "MaximalAbelianGaugeFixing" << std::endl
			<<  "MaximalAbelianProjection" << std::endl
			<<  "LandauGluonPropagator" << std::endl
			<<  "LandauGhostPropagator" << std::endl
			<<  "ScalarFermionHMC" << std::endl
			<<  "HiggsGaugeHMC" << std::endl
			<<  "RandomScalarInitializer" << std::endl
			<<  "MCScalar" << std::endl
			<<  "RandomGaugeTransformation" << std::endl
			<<  "MeanScalarField" << std::endl;
		}

	}

} /* namespace Update */
