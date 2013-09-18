/*
 * LatticeSweep.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: spiem_01
 */

#include "LatticeSweep.h"
#include "PureGaugeUpdater.h"
#include "Plaquette.h"
#include "FromString.h"
#include "PureGaugeHMCUpdater.h"
#include "TwoFlavorHMCUpdater.h"
#include "NFlavorQCDUpdater.h"
#include "NFlavorBlockUpdater.h"
#include "PureGaugeOverrelaxation.h"
#include "OutputSweep.h"
#include "TestLinearAlgebra.h"
#include "Eigenvalues.h"
#include "MesonCorrelator.h"
#include "ChiralCondensate.h"
#include "PolyakovLoop.h"
#include "Glueball.h"
#include "ReUnit.h"
#include "BandTwoFlavorUpdater.h"
#include "PureGaugeWilsonLoops.h"
#include "TestCommunication.h"
#include "GluinoGlue.h"

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
	} else if (name == "Eigenvalues") {
		return new Eigenvalues();
	} else if (name == "MesonCorrelator") {
		return new MesonCorrelator();
	} else if (name == "ChiralCondensate") {
		return new ChiralCondensate();
	} else if (name == "PolyakovLoop") {
		return new PolyakovLoop();
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
	} else if (name == "TestCommunication") {
		return new TestCommunication();
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

void LatticeSweep::printSweepsName() {
	if (isOutputProcess()) {
		std::cout << "List of the possible Sweeps name:" << std::endl;
		std::cout
		<< " PureGaugeCM" << std::endl
		<< " Plaquette" << std::endl
		<< " PureGaugeHMC" << std::endl
		<< " TwoFlavorQCD" << std::endl
		<< " NFlavorQCD" << std::endl
		<< " NFlavorBlock" << std::endl
		<< " PureGaugeOverrelaxation" << std::endl
		<< " Output" << std::endl
		<< " TestLinearAlgebra" << std::endl
		<< " Eigenvalues" << std::endl
		<< " MesonCorrelator" << std::endl
		<< " ChiralCondensate" << std::endl
		<< " PolyakovLoop" << std::endl
		<< " Glueball" << std::endl
		<< " ReUnit" << std::endl
		<< " BandTwoFlavorHMCUpdater" << std::endl
		<< " PureGaugeWilsonLoops" << std::endl
		<< " TestCommunication" << std::endl;
	}

}

} /* namespace Update */
