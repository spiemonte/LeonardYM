/*
 * ComplementBlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "SAPPreconditioner.h"
#include "algebra_utils/AlgebraUtils.h"

namespace Update {

SAPPreconditioner::SAPPreconditioner(DiracOperator* _diracOperator, ComplementBlockDiracOperator* _K) : DiracOperator(), diracOperator(_diracOperator), K(_K), steps(7), precision(0.00001) { }

SAPPreconditioner::~SAPPreconditioner() { }

void SAPPreconditioner::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	//Adaptative strategy
	long_real_t norm = AlgebraUtils::squaredNorm(input);
	if (norm < 10.*precision) K->setPrecision(0.001*norm);
	else K->setPrecision(precision);

	output = input;
	K->multiply(tmp1,input);
	for (int i = 0; i < steps; ++i) {
		diracOperator->multiply(tmp2,output);
		//TODO Use the same vector?
		K->multiply(tmp3,tmp2);
#pragma omp parallel for
		for (int site = 0; site < tmp2.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] = output[site][mu] - tmp3[site][mu] + tmp1[site][mu];
			}
		}
	}
}

void SAPPreconditioner::multiplyAdd(reduced_dirac_vector_t& , const reduced_dirac_vector_t& , const reduced_dirac_vector_t& , const std::complex<real_t>& ) {
	//TODO: to be implemented, not needed for computation of propagators
	
}

void SAPPreconditioner::setLattice(const extended_fermion_lattice_t& _lattice) {
	this->lattice = _lattice;
	diracOperator->setLattice(_lattice);
	K->setLattice(_lattice);
}

FermionForce* SAPPreconditioner::getForce() const {
	return diracOperator->getForce();
}

void SAPPreconditioner::setKappa(real_t _kappa) {
	kappa = _kappa;
	diracOperator->setKappa(_kappa);
	K->setKappa(_kappa);
}

void SAPPreconditioner::setSteps(int _steps) {
	steps = _steps;
}

void SAPPreconditioner::setPrecision(double _precision) {
	precision = _precision;
	K->setPrecision(_precision);
}

} /* namespace Update */
