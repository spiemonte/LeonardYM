/*
 * DiracOperator.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#include "DiracOperator.h"
#include "OverlapOperator.h"
#include "ExactOverlapOperator.h"
#include "SquareOverlapOperator.h"
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

DiracOperator* DiracOperator::getInstance(const std::string& name, unsigned int power, const StorageParameters& parameters, const std::string& basename) {
	if (power == 1) {
		if (name == "DiracWilson") {
			DiracWilsonOperator* result = new DiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->name = name;
			return result;
		}
		else if (name == "BasicDiracWilson") {
			BasicDiracWilsonOperator* result = new BasicDiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->name = name;
			return result;
		}
		else if (name == "Improved") {
			ImprovedDiracWilsonOperator* result = new ImprovedDiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setCSW(parameters.get<double>(basename+"csw"));
			result->name = name;
			return result;
		}
		else if (name == "Overlap") {
			OverlapOperator* result = new OverlapOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setMass(parameters.get<double>(basename+"mass"));

			std::vector< std::complex<real_t> > coeff = parameters.get< std::vector< std::complex<real_t> > >(basename+"OverlapOperator::squareRootApproximation");
			Polynomial squareRootApproximation;
			squareRootApproximation.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			squareRootApproximation.setRoots(coeff);


			result->setSquareRootApproximation(squareRootApproximation);
			result->name = name;
			return result;
		}
		else if (name == "ExactOverlap") {
			ExactOverlapOperator* result = new ExactOverlapOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setMass(parameters.get<double>(basename+"mass"));

			std::vector< std::complex<real_t> > coeff = parameters.get< std::vector< std::complex<real_t> > >(basename+"ExactOverlapOperator::squareRootApproximation");
			Polynomial squareRootApproximation;
			squareRootApproximation.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			squareRootApproximation.setRoots(coeff);

			result->setSquareRootApproximation(squareRootApproximation);
			result->name = name;

			result->getDiracEigenSolver()->setTolerance(parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::eigensolver_precision"));
			if (parameters.get<std::string>(basename+"ExactOverlapOperator::eigensolver::use_chebyshev") == "true") {
				unsigned int chebyshevOrder = parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::chebyshev_order");
				double chebyshevLeft = parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::chebyshev_left");
				double chebyshevRight = parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::chebyshev_right");
				result->getDiracEigenSolver()->setUseChebyshev(true);
				result->getDiracEigenSolver()->getChebyshevRecursion()->setParameters(chebyshevLeft, chebyshevRight, chebyshevOrder);
			}

			result->getDiracEigenSolver()->setMaximalNumberOfRestarts(parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver"));
			result->setNumberOfEigenvalues(parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::number_eigenvalues"));
			return result;
		}
		else {
			std::cout << "Dirac Wilson Operator" << name << " not supported!" << std::endl;
			exit(1);
		}
	} else if (power == 2) {
		if (name == "DiracWilson") {
			SquareDiracWilsonOperator* result = new SquareDiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->name = name;
			return result;
		}
		else if (name == "Improved") {
			SquareImprovedDiracWilsonOperator* result = new SquareImprovedDiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setCSW(parameters.get<double>(basename+"csw"));
			result->name = name;
			return result;
		}
		else if (name == "BasicDiracWilson") {
			BasicSquareDiracWilsonOperator* result = new BasicSquareDiracWilsonOperator();
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->name = name;
			return result;
		}
		else if (name == "Overlap") {
			OverlapOperator* ov = new OverlapOperator();
			SquareOverlapOperator* result = new SquareOverlapOperator();
			result->setOverlapOperator(ov);
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setMass(parameters.get<double>(basename+"mass"));

			std::vector< std::complex<real_t> > coeff = parameters.get< std::vector< std::complex<real_t> > >(basename+"OverlapOperator::squareRootApproximation");
			Polynomial squareRootApproximation;
			squareRootApproximation.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			squareRootApproximation.setRoots(coeff);


			result->setSquareRootApproximation(squareRootApproximation);
			result->name = name;
			return result;
		}
		else if (name == "ExactOverlap") {
			ExactOverlapOperator* ov = new ExactOverlapOperator();
			
			ov->getDiracEigenSolver()->setTolerance(parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::eigensolver_precision"));
			if (parameters.get<std::string>(basename+"ExactOverlapOperator::eigensolver::use_chebyshev") == "true") {
				unsigned int chebyshevOrder = parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::chebyshev_order");
				double chebyshevLeft = parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::chebyshev_left");
				double chebyshevRight = parameters.get<double>(basename+"ExactOverlapOperator::eigensolver::chebyshev_right");
				ov->getDiracEigenSolver()->setUseChebyshev(true);
				ov->getDiracEigenSolver()->getChebyshevRecursion()->setParameters(chebyshevLeft, chebyshevRight, chebyshevOrder);
			}

			ov->getDiracEigenSolver()->setMaximalNumberOfRestarts(parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver"));
			ov->setNumberOfEigenvalues(parameters.get<unsigned int>(basename+"ExactOverlapOperator::eigensolver::number_eigenvalues"));

			SquareOverlapOperator* result = new SquareOverlapOperator();
			result->setOverlapOperator(ov);
			result->setKappa(parameters.get<double>(basename+"kappa"));
			result->setMass(parameters.get<double>(basename+"mass"));

			std::vector< std::complex<real_t> > coeff = parameters.get< std::vector< std::complex<real_t> > >(basename+"ExactOverlapOperator::squareRootApproximation");
			Polynomial squareRootApproximation;
			squareRootApproximation.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			squareRootApproximation.setRoots(coeff);


			result->setSquareRootApproximation(squareRootApproximation);
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

void DiracOperator::registerParameters(po::options_description& desc, const std::string& basename) {
	desc.add_options()
		((basename+"dirac_operator").c_str(), po::value<std::string>()->default_value("DiracWilson"), "The name of the Dirac operator")
		((basename+"kappa").c_str(), po::value<double>()->default_value(0.0), "set the value of the kappa for the Dirac Wilson operator")
		((basename+"csw").c_str(), po::value<double>()->default_value(0.0), "set the value of the clover term for the Dirac Wilson operator")
		((basename+"mass").c_str(), po::value<double>()->default_value(0.0), "set the value of the mass for the Overlap operator")
		((basename+"OverlapOperator::squareRootApproximation").c_str(), po::value<std::string>()->default_value(""), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		((basename+"ExactOverlapOperator::squareRootApproximation").c_str(), po::value<std::string>()->default_value(""), "Approximation of x^(1/2) used for the Overlap fermion sign function (syntax: {(scalingre,scalingim),(r1re,r1im), ..., (rnre,rnim)})")
		((basename+"ExactOverlapOperator::eigensolver::use_chebyshev").c_str(), po::value<std::string>()->default_value("false"), "Use Chebyshev acceleration? (true/false)")
		((basename+"ExactOverlapOperator::eigensolver::chebyshev_left").c_str(), po::value<double>()->default_value(0.2), "Left interval of the Chebyshev polynomial.")
		((basename+"ExactOverlapOperator::eigensolver::chebyshev_right").c_str(), po::value<double>()->default_value(7.), "Right interval of the Chebyshev polynomial.")
		((basename+"ExactOverlapOperator::eigensolver::chebyshev_order").c_str(), po::value<unsigned int>()->default_value(15), "Order of the Chebyshev acceleration. It must be an odd number")
		((basename+"ExactOverlapOperator::eigensolver::eigensolver_precision").c_str(), po::value<double>()->default_value(0.000000001), "set the precision used by the eigensolver")
		((basename+"ExactOverlapOperator::eigensolver::number_extra_vectors").c_str(), po::value<unsigned int>(), "Number of extra vectors for the Arnoldi algorithm used in the computation of the eigenvectors, increase this number to increase precision")
		((basename+"ExactOverlapOperator::eigensolver::maximal_number_restarts_eigensolver").c_str(), po::value<unsigned int>()->default_value(50), "Number of restarts for the implicitly restarted Arnoldi algorithm")
		((basename+"ExactOverlapOperator::eigensolver::number_eigenvalues").c_str(), po::value<unsigned int>(), "Number of eigenvalues of the dirac wilson operator to be computed")
		;
}

DiracOperator* DiracOperator::getSquareRoot(DiracOperator* dirac) {
	if (dirac->name == "DiracWilson") {
		if (dynamic_cast<SquareDiracWilsonOperator*>(dirac)) {
			DiracWilsonOperator* result = new DiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->lattice = dirac->lattice;
			result->gamma5 = dirac->gamma5;
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
			result->gamma5 = dirac->gamma5;
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
			result->gamma5 = dirac->gamma5;
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

DiracOperator* DiracOperator::getSquare(DiracOperator* dirac) {
	if (dirac->name == "DiracWilson") {
		if (dynamic_cast<DiracWilsonOperator*>(dirac)) {
			SquareDiracWilsonOperator* result = new SquareDiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->setLattice(dirac->lattice);
			result->gamma5 = dirac->gamma5;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "Improved") {
		if (dynamic_cast<ImprovedDiracWilsonOperator*>(dirac)) {
			SquareImprovedDiracWilsonOperator* result = new SquareImprovedDiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->setCSW(dynamic_cast<ImprovedDiracWilsonOperator*>(dirac)->getCSW());
			result->name = dirac->name;
			result->setLattice(dirac->lattice);
			result->gamma5 = dirac->gamma5;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "BasicDiracWilson") {
		if (dynamic_cast<BasicDiracWilsonOperator*>(dirac)) {
			BasicSquareDiracWilsonOperator* result = new BasicSquareDiracWilsonOperator();
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->setLattice(dirac->lattice);
			result->gamma5 = dirac->gamma5;
			return result;
		} else {
			std::cout << "Power of the Dirac Wilson Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "Overlap") {
		if (dynamic_cast<OverlapOperator*>(dirac)) {
			OverlapOperator* op = dynamic_cast<OverlapOperator*>(dirac);
			SquareOverlapOperator* result = new SquareOverlapOperator();
			result->setOverlapOperator(op);
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->setLattice(dirac->lattice);
			result->gamma5 = dirac->gamma5;
			result->setMass(op->getMass());
			result->setSquareRootApproximation(op->getSquareRootApproximation());
			return result;
		} else {
			std::cout << "Power of the Overlap Operator not supported!" << std::endl;
			exit(1);
		}
	}
	else if (dirac->name == "ExactOverlap") {
		if (dynamic_cast<OverlapOperator*>(dirac)) {
			ExactOverlapOperator* op = dynamic_cast<ExactOverlapOperator*>(dirac);
			SquareOverlapOperator* result = new SquareOverlapOperator();
			result->setOverlapOperator(new ExactOverlapOperator(*op));
			result->setKappa(dirac->getKappa());
			result->name = dirac->name;
			result->setLattice(dirac->lattice);
			result->gamma5 = dirac->gamma5;
			result->setMass(op->getMass());
			result->setSquareRootApproximation(op->getSquareRootApproximation());
			return result;
		} else {
			std::cout << "Power of the Overlap Operator not supported!" << std::endl;
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
