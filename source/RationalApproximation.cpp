/*
 * RationalApproximation.cpp
 *
 *  Created on: May 10, 2012
 *      Author: spiem_01
 */

#include "RationalApproximation.h"
#include "MEMultishiftSolver.h"
#include "MMMRMultishiftSolver.h"
#include "AlgebraUtils.h"

namespace Update {

std::vector< extended_dirac_vector_t > RationalApproximation::tmp;
//extended_dirac_vector_t RationalApproximation::tmp1;
//extended_dirac_vector_t RationalApproximation::tmp2;

RationalApproximation::RationalApproximation() : precision(0.0000000000001), preconditionerPrecision(0.00000000000001), maximumSteps(3000), multishiftSolver(new MMMRMultishiftSolver()), biConjugateGradient(new BiConjugateGradient()), preconditionerRecursion(0) { }

RationalApproximation::RationalApproximation(const RationalApproximation& toCopy) : alphas(toCopy.alphas), betas(toCopy.betas), precision(toCopy.precision),  preconditionerPrecision(toCopy.preconditionerPrecision), maximumSteps(toCopy.maximumSteps), multishiftSolver(new MMMRMultishiftSolver()), biConjugateGradient(new BiConjugateGradient()), preconditionerRecursion(toCopy.preconditionerRecursion) { }

RationalApproximation::RationalApproximation(const std::vector< real_t >& _alphas, const std::vector< real_t >& _betas) : alphas(_alphas), betas(_betas), precision(0.0000000000001), preconditionerPrecision(0.00000000000001), maximumSteps(3000), multishiftSolver(new MMMRMultishiftSolver()), biConjugateGradient(new BiConjugateGradient()), preconditionerRecursion(0) { }

RationalApproximation::~RationalApproximation() {
	delete multishiftSolver;
	delete biConjugateGradient;
}

void RationalApproximation::evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& output, const extended_dirac_vector_t& input, DiracOperator* preconditioner) {
	//Allocate the memory
	if (tmp.size() != betas.size()) tmp.resize(betas.size());
	//get and set the multishiftSolver
	multishiftSolver->setPrecision(precision);
	multishiftSolver->setMaxSteps(maximumSteps);
	biConjugateGradient->setMaximumSteps(maximumSteps);
	biConjugateGradient->setPrecision(preconditionerPrecision);
	tmp1 = input;
	if (preconditioner != 0) {
		//First we apply the transformation
		for (unsigned int k = 0; k < preconditionerRecursion; ++k) {
			preconditioner->multiply(output,tmp1);
			tmp1 = output;
		}
	}
	//Solve the dirac equations
	multishiftSolver->solve(diracOperator, tmp1, tmp, betas);
	//Test of multishift
	/*std::vector<real_t>::iterator beta;
	std::vector<dirac_vector_t>::iterator vv;
	dirac_vector_t ll;
	for (beta = betas.begin(), vv = tmp.begin(); beta != betas.end(); ++beta, ++vv) {
		diracOperator->multiplyAdd(ll,*vv,*vv,*beta);
		if (isOutputProcess()) std::cout << "Test for beta" << *beta << std::endl;
		if (isOutputProcess()) std::cout << ll[5][3] << std::endl;
		if (isOutputProcess()) std::cout << input[5][3] << std::endl;
	}*/
	//Compute the result
	std::vector<real_t>::iterator alpha;
	std::vector<extended_dirac_vector_t>::iterator vector;
	for (int site = 0; site < output.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(output[site][mu]);
		}
		for (alpha = alphas.begin(), vector = tmp.begin(); vector != tmp.end(); ++vector, ++alpha) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] += (*alpha)*(*vector)[site][mu];
			}
		}
	}
	output.updateHalo();
	if (preconditioner != 0) {
		unsigned int invertSteps = 0;
		//Now we invert the initial transformation
		preconditioner->setGamma5(false);
		for (unsigned int k = 0; k < preconditionerRecursion; ++k) {
			//We apply first gamma_5
#pragma omp parallel for
			for (int site = 0; site < output.completesize; ++site) {
				output[site][2] = -output[site][2];
				output[site][3] = -output[site][3];
			}
			biConjugateGradient->solve(preconditioner, output, tmp1);
			invertSteps += biConjugateGradient->getLastSteps();
			output = tmp1;
		}
		preconditioner->setGamma5(true);
		if (preconditionerRecursion != 0 && isOutputProcess()) std::cout << "RationalApproximation::Map convergence in " << invertSteps << " steps." << std::endl;
	}
}

complex RationalApproximation::evaluate(const complex& x) const {
	std::vector<real_t>::const_iterator alpha;
	std::vector<real_t>::const_iterator beta;
	complex result = 0.;
	for (alpha = alphas.begin(), beta = betas.begin(); alpha != alphas.end(); ++alpha, ++beta) {
		result += (*alpha)/(x + (*beta));
	}
	return result;	
}

void RationalApproximation::setAlphas(const std::vector< real_t >& _alphas) {
	alphas = _alphas;
}

void RationalApproximation::setBetas(const std::vector< real_t >& _betas) {
	betas = _betas;
}

std::vector< real_t > RationalApproximation::getAlphas() const {
	return alphas;
}

std::vector< real_t > RationalApproximation::getBetas() const {
	return betas;
}

void RationalApproximation::setPrecision(const real_t& _precision) {
	precision = _precision;
}

real_t RationalApproximation::getPrecision() const {
	return precision;
}

void RationalApproximation::setPreconditionerPrecision(const real_t& _precision) {
	preconditionerPrecision = _precision;
}

real_t RationalApproximation::getPreconditionerPrecision() const {
	return preconditionerPrecision;
}

void RationalApproximation::setMaximumRecursion(unsigned int maximumRecursion) {
	maximumSteps = maximumRecursion;
}

unsigned int RationalApproximation::getMaximumRecursion() const {
	return maximumSteps;
}

void RationalApproximation::setPreconditionerRecursion(unsigned int recursion) {
	preconditionerRecursion = recursion;
}

unsigned int RationalApproximation::getPreconditionerRecursion() const {
	return preconditionerRecursion;
}

} /* namespace Update */
