/*
 * MultishiftSolver.cpp
 *
 *  Created on: May 11, 2012
 *      Author: spiem_01
 */

#include "MultishiftSolver.h"
#include "MEMultishiftSolver.h"
#include "MMMRMultishiftSolver.h"
#include "ChronologicalMultishiftSolver.h"
#define MULTISHIFTLOG

namespace Update {

MultishiftSolver::MultishiftSolver(real_t _epsilon, unsigned int _maxSteps) : epsilon(_epsilon), maxSteps(_maxSteps) { }

MultishiftSolver::~MultishiftSolver() { }

MultishiftSolver* MultishiftSolver::getInstance(const std::string& name) {
	if (name == "mass_estrapolation") {
		return new MEMultishiftSolver();
	}
	else if (name == "chronological_mass_estrapolation") {
		return new ChronologicalMultishiftSolver();
	}
	else if (name == "minimal_residual") {
		return new MMMRMultishiftSolver();
	}
	else {
		if (isOutputProcess()) std::cout << "Name " << name << " of multishift solver is not recognized!" << std::endl;
		exit(1);
	}
}


void MultishiftSolver::setPrecision(double _epsilon) {
	epsilon = _epsilon;
}

double MultishiftSolver::getPrecision() const {
	return epsilon;
}

void MultishiftSolver::setMaxSteps(unsigned int _maxSteps) {
	maxSteps = _maxSteps;
}

unsigned int MultishiftSolver::getMaxSteps() const {
	return maxSteps;
}

} /* namespace Update */
