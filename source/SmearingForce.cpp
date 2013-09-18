/*
 * SmearingForce.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: spiem_01
 */

#include "SmearingForce.h"

namespace Update {

SmearingForce::SmearingForce() { }

SmearingForce::~SmearingForce() { }

FermionicForceMatrix SmearingForce::derivative(environment_t& env, const FermionicForceMatrix& unsmearedDerivative, int site, unsigned int mu) {
	for (unsigned int level = env.configurations.get<unsigned int>("level_smearing_updater"); level > 0; --level) {
		env.smearedConfigurations[level];
	}


	adjoint_lattice_t& prevLattice = env.smearedConfigurations[0];
	//This tensor contains the derivative \frac{\partial U_{ij}}{\partial U_{kl}}
	std::complex<real_t> tensorDerivative[diracVectorLength][diracVectorLength][diracVectorLength][diracVectorLength];
	//Set to zero
	for (unsigned int i = 0; i < diracVectorLength; ++i) {
		for (unsigned int j = 0; j < diracVectorLength; ++j) {
			for (unsigned int k = 0; k < diracVectorLength; ++k) {
				for (unsigned int l = 0; l < diracVectorLength; ++l) {
					tensorDerivative[i][j][k][l] = 0.;
				}
			}
		}
	}

	//The first term that we add is the derivative of the link off the exponential
	WilsonGaugeAction action;
	FermionicGroup staple = ;
	FermionicGroup Q = ;
	FermionicGroup exponential = ;
	for (unsigned int i = 0; i < diracVectorLength; ++i) {
		for (unsigned int j = 0; j < diracVectorLength; ++j) {
			for (unsigned int k = 0; k < diracVectorLength; ++k) {
				for (unsigned int l = 0; l < diracVectorLength; ++l) {
					if (j == l) tensorDerivative[i][j][k][l] += exponential.at(i,k);
				}
			}
		}
	}

	//Then we add the derivative of the exponential
	//For doing this, we need the derivative of the single coefficients
	for (unsigned int i = 0; i < diracVectorLength; ++i) {
		for (unsigned int j = 0; j < diracVectorLength; ++j) {
			for (unsigned int k = 0; k < diracVectorLength; ++k) {
				for (unsigned int l = 0; l < diracVectorLength; ++l) {
					if (j == l) tensorDerivative[i][j][k][l] += exponential.at(i,k);
				}
			}
		}
	}

	std::complex<real_t> q0 = Q(0,1);
	std::complex<real_t> q1 = Q(0,2);
	std::complex<real_t> q2 = Q(1,2);

	std::complex<real_t> df1_q[diracVectorLength];
	std::complex<real_t> df2_q[diracVectorLength];

	//The su2 norm
	std::complex<real_t> norm = sqrt(-q0*q0 - q1*q1 - q2*q2);
	//the prefactor common the first derivatives
	std::complex<real_t> prefactor_f1 = - exp(-norm)/(2.*norm*norm*norm) + exp(norm)/(2.*norm*norm*norm) - exp(-norm)/(2.*norm*norm) + exp(norm)/(2.*norm*norm);

	//This is the derivative of f1 respect to q_i
	df1_q[0] = q0*prefactor_f1;
	df1_q[1] = q1*prefactor_f1;
	df1_q[2] = q2*prefactor_f1;

	//the prefactor common the first derivatives
	std::complex<real_t> prefactor_f2 = -(2.)/(norm*norm) + exp(-norm)/(norm*norm) + exp(norm)/(norm*norm) + exp(-norm)/(2.*norm*norm*norm) - exp(norm)/(2.*norm*norm*norm);

	//This is the derivative of f1 respect to q_i
	df2_q[0] = q0*prefactor_f2;
	df2_q[1] = q1*prefactor_f2;
	df2_q[2] = q2*prefactor_f2;

	std::complex<real_t> dexp[diracVectorLength][diracVectorLength][diracVectorLength][diracVectorLength];

	//This terms contains the derivative of q respect to the staple
	//Q = S - Transpose(S)
	std::complex<real_t> dq0[diracVectorLength][diracVectorLength];
	std::complex<real_t> dq1[diracVectorLength][diracVectorLength];
	std::complex<real_t> dq2[diracVectorLength][diracVectorLength];

	for (unsigned int i = 0; i < diracVectorLength; ++i) {
		for (unsigned int j = 0; j < diracVectorLength; ++j) {
			if (i == 1 && j == 2) dq0[i][j] = 1.;
			else if (i == 2 && j == 1) dq0[i][j] = -1.;
			else dq0[i][j] = 0.;

			if (i == 1 && j == 3) dq1[i][j] = 1.;
			else if (i == 3 && j == 1) dq1[i][j] = -1.;
			else dq1[i][j] = 0.;

			if (i == 2 && j == 3) dq0[i][j] = 1.;
			else if (i == 3 && j == 2) dq0[i][j] = -1.;
			else dq2[i][j] = 0.;
		}
	}



	for (unsigned int i = 0; i < diracVectorLength; ++i) {
		for (unsigned int j = 0; j < diracVectorLength; ++j) {

		}
	}

}



} /* namespace Update */
