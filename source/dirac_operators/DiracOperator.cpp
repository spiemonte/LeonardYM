/*
 * DiracOperator.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#include "DiracOperator.h"
#include "DiracWilsonOperator.h"
#include "ImprovedDiracWilsonOperator.h"
#include "SquareDiracWilsonOperator.h"
#include "SquareImprovedDiracWilsonOperator.h"
#include "BasicDiracWilsonOperator.h"
#include "BasicSquareDiracWilsonOperator.h"

namespace Update {

DiracOperator::DiracOperator() : kappa(0.), gamma5(true) { }

DiracOperator::DiracOperator(const extended_fermion_lattice_t& _lattice, double _kappa, bool _gamma5) : lattice(_lattice), kappa(_kappa), gamma5(_gamma5) { }

DiracOperator::~DiracOperator() { }

DiracOperator* DiracOperator::getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters) {
	if (power == 1) {
		if (name == "DiracWilson") {
			DiracWilsonOperator* result = new DiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->name = name;
			return result;
		}
		else if (name == "BasicDiracWilson") {
			BasicDiracWilsonOperator* result = new BasicDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->name = name;
			return result;
		}
		else if (name == "Improved") {
			ImprovedDiracWilsonOperator* result = new ImprovedDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->setCSW(parameters.get<double>("csw"));
			result->name = name;
			return result;
		}
		else {
			std::cout << "Dirac Wilson Operator" << name << " not supported!" << std::endl;
			exit(1);
		}
	} else if (power == 2) {
		if (name == "DiracWilson") {
			SquareDiracWilsonOperator* result = new SquareDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->name = name;
			return result;
		}
		else if (name == "Improved") {
			SquareImprovedDiracWilsonOperator* result = new SquareImprovedDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->setCSW(parameters.get<double>("csw"));
			result->name = name;
			return result;
		}
		else if (name == "BasicDiracWilson") {
			BasicSquareDiracWilsonOperator* result = new BasicSquareDiracWilsonOperator();
			result->setKappa(parameters.get<double>("kappa"));
			result->name = name;
			return result;
		}
		else {
			std::cout << "Dirac Wilson Operator" << name << " not supported!" << std::endl;
			exit(1);
		}
	} else {
		std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
		exit(1);
	}
}

DiracOperator* DiracOperator::getSquareRoot(DiracOperator* dirac) {
	if (dirac->name == "DiracWilson") {
		if (dynamic_cast<SquareDiracWilsonOperator*>(dirac)) {
			DiracWilsonOperator* result = new DiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->lattice = dirac->lattice;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "Improved") {
		if (dynamic_cast<SquareImprovedDiracWilsonOperator*>(dirac)) {
			ImprovedDiracWilsonOperator* result = new ImprovedDiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->setCSW(dynamic_cast<SquareImprovedDiracWilsonOperator*>(dirac)->getCSW());
			result->name = dirac->name;
			result->lattice = dirac->lattice;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "BasicDiracWilson") {
		if (dynamic_cast<BasicSquareDiracWilsonOperator*>(dirac)) {
			BasicDiracWilsonOperator* result = new BasicDiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->lattice = dirac->lattice;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else {
		std::cout << "Dirac Wilson Operator" << dirac->name << " not supported!" << std::endl;
		exit(1);
	}
}

void DiracOperator::setKappa(real_t _kappa) {
	kappa = _kappa;
}

real_t DiracOperator::getKappa() const {
	return kappa;
}

void DiracOperator::setGamma5(bool _gamma5) {
	gamma5 = _gamma5;
}

bool DiracOperator::getGamma5() const {
	return gamma5;
}

const reduced_fermion_lattice_t *DiracOperator::getLattice() const {
	return &lattice;
}

void DiracOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
}

std::string DiracOperator::getName() const {
	return name;
}

} /* namespace Update */
