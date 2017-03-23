/*
 * LatticeSweep.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#include "LatticeSweep.h"
#include "pure_gauge/PureGaugeUpdater.h"
#include "wilson_loops/Plaquette.h"
#include "utils/FromString.h"
#include "hmc_updaters/PureGaugeHMCUpdater.h"
#include "hmc_updaters/TwoFlavorHMCUpdater.h"
#include "hmc_updaters/NFlavorQCDUpdater.h"
#include "hmc_updaters/NFlavorBlockUpdater.h"
#include "pure_gauge/PureGaugeOverrelaxation.h"
#include "io/OutputSweep.h"
#include "tests/TestLinearAlgebra.h"
#include "tests/TestSpeedDiracOperators.h"
#include "fermion_measurements/Eigenvalues.h"
#include "correlators/MesonCorrelator.h"
#include "fermion_measurements/ChiralCondensate.h"
#include "fermion_measurements/SingletOperators.h"
#include "fermion_measurements/XSpaceCorrelators.h"
#include "fermion_measurements/NPRVertex.h"
#include "polyakov_loops/PolyakovLoop.h"
#include "correlators/Glueball.h"
#include "utils/ReUnit.h"
#include "utils/LandauGaugeFixing.h"
#include "utils/MaximalAbelianGaugeFixing.h"
#include "utils/MaximalAbelianProjection.h"
#include "utils/LandauPropagator.h"
#include "hmc_updaters/BandTwoFlavorUpdater.h"
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
#include "hmc_updaters/MultiStepNFlavorQCDUpdater.h"
#include "hmc_updaters/TwistedMultiStepNFlavorQCDUpdater.h"

namespace Update {

LatticeSweep::LatticeSweep() { }

LatticeSweep::LatticeSweep(unsigned int _numberTimes, unsigned int _sweepToJump) : numberTimes(_numberTimes), sweepToJump(_sweepToJump) { }

LatticeSweep::~LatticeSweep() { }

LatticeSweep* LatticeSweep::getInstance(const std::string& name) {
	if (name == "PureGaugeCM") {
		return new PureGaugeUpdater();
	} else if (name == "Plaquette") {
		return new Plaquette();
	} else if (name == "PureGaugeHMC") {
		return new PureGaugeHMCUpdater();
	} else if (name == "TwoFlavorQCD") {
		return new TwoFlavorHMCUpdater();
	} else if (name == "NFlavorQCD") {
		return new NFlavorQCDUpdater();
	} else if (name == "NFlavorBlock") {
		return new NFlavorBlockUpdater();
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
	} else if (name == "BandTwoFlavorHMCUpdater") {
		return new BandTwoFlavorUpdater();
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
	} else if (name == "MultiStepNFlavorQCD") {
		return new MultiStepNFlavorQCDUpdater();
	} else if (name == "TwistedMultiStepNFlavorQCD") {
		return new TwistedMultiStepNFlavorQCDUpdater();
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
	} else if (name == "LandauPropagator") {
		return new LandauPropagator();
	}
	else {
		if (isOutputProcess()) std::cout << "Unknown name sweep: " << name << std::endl;
		exit(1);
	}
}

void LatticeSweep::call(environment_t& environment) {
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

void LatticeSweep::addParameters(po::options_description& desc) {
	PureGaugeUpdater::registerParameters(desc);
	Plaquette::registerParameters(desc);
	PureGaugeHMCUpdater::registerParameters(desc);
	TwoFlavorHMCUpdater::registerParameters(desc);
	NFlavorQCDUpdater::registerParameters(desc);
	NFlavorBlockUpdater::registerParameters(desc);
	PureGaugeOverrelaxation::registerParameters(desc);
	OutputSweep::registerParameters(desc);
	TestLinearAlgebra::registerParameters(desc);
	TestSpeedDiracOperators::registerParameters(desc);
	Eigenvalues::registerParameters(desc);
	MesonCorrelator::registerParameters(desc);
	ChiralCondensate::registerParameters(desc);
	PolyakovLoop::registerParameters(desc);
	PolyakovLoopEigenvalues::registerParameters(desc);
	PolyakovLoopCorrelator::registerParameters(desc);
	AdjointPolyakovLoop::registerParameters(desc);
	Glueball::registerParameters(desc);
	GluinoGlue::registerParameters(desc);
	ReUnit::registerParameters(desc);
	BandTwoFlavorUpdater::registerParameters(desc);
	PureGaugeWilsonLoops::registerParameters(desc);
	WilsonLoop::registerParameters(desc);
	TestCommunication::registerParameters(desc);
	ReadGaugeConfiguration::registerParameters(desc);
	WilsonFlow::registerParameters(desc);
	GaugeEnergy::registerParameters(desc);
	MultiStepNFlavorQCDUpdater::registerParameters(desc);
	TwistedMultiStepNFlavorQCDUpdater::registerParameters(desc);
	SingletOperators::registerParameters(desc);
	XSpaceCorrelators::registerParameters(desc);
	NPRVertex::registerParameters(desc);
	LandauGaugeFixing::registerParameters(desc);
	MaximalAbelianGaugeFixing::registerParameters(desc);
	MaximalAbelianProjection::registerParameters(desc);
	LandauPropagator::registerParameters(desc);
}

void LatticeSweep::registerParameters(po::options_description&) {
}

void LatticeSweep::printSweepsName() {
	if (isOutputProcess()) {
		std::cout << "List of the possible Sweeps name:" << std::endl;
		std::cout
		<< "PureGaugeCM" << std::endl
		<< "Plaquette" << std::endl
		<< "PureGaugeHMC" << std::endl
		<< "TwoFlavorQCD" << std::endl
		<< "NFlavorQCD" << std::endl
		<< "NFlavorBlock" << std::endl
		<< "PureGaugeOverrelaxation" << std::endl
		<< "Output" << std::endl
		<< "TestLinearAlgebra" << std::endl
		<< "TestSpeedDiracOperators" << std::endl
		<< "Eigenvalues" << std::endl
		<< "MesonCorrelator" << std::endl
		<< "ChiralCondensate" << std::endl
		<< "PolyakovLoop" << std::endl
		<< "AdjointPolyakovLoop" << std::endl
		<< "Glueball" << std::endl
		<< "GluinoGlue" << std::endl
		<< "ReUnit" << std::endl
		<< "BandTwoFlavorHMCUpdater" << std::endl
		<< "PureGaugeWilsonLoops" << std::endl
		<< "WilsonLoop" << std::endl
		<< "TestCommunication" << std::endl
		<< "ReadGaugeConfiguration" << std::endl
		<< "WilsonFlow" << std::endl
		<< "GaugeEnergy" << std::endl
		<< "MultiStepNFlavorQCD" << std::endl
		<< "TwistedMultiStepNFlavorQCD" << std::endl
		<< "LandauPropagator" << std::endl;
	}

}

} /* namespace Update */
