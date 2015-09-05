/*
 * SmearingForce.cpp
 *
 *  Created on: Dec 7, 2012
 *      Author: spiem_01
 */

#include "SmearingForce.h"
#define REDUCTION 1.4
#define MAXSIZE 25

namespace Update {

real_t fabs(const std::complex<real_t>& number) {
	return real(conj(number)*number);
}

SmearingForce::SmearingForce() : StoutSmearing() { }

SmearingForce::~SmearingForce() { }

void SmearingForce::derivative(const extended_force_lattice_t& smearedDerivative, extended_force_lattice_t& unsmearedDerivative, extended_gauge_lattice_t& unsmearedLattice, real_t rho) {
	typedef extended_gauge_lattice_t LT;
	
#pragma omp parallel for
	for (int site = 0; site < unsmearedLattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			//First term
			for (int color = 0; color < numberColors*numberColors - 1; ++color) {
				unsmearedDerivative[site][mu][color] = smearedDerivative[site][mu]*this->ridder(unsmearedLattice, site, mu, color, site, mu, rho);
			}

			//Second term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[LT::sup(site,mu)][nu]*this->ridder(unsmearedLattice, site, mu, color, LT::sup(site,mu), nu, rho);
					}
				}
			}

			//Third term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[LT::sup(site,nu)][mu]*this->ridder(unsmearedLattice, site, mu, color, LT::sup(site,nu), mu, rho);
					}
				}
			}

			//Fourth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[site][nu]*this->ridder(unsmearedLattice, site, mu, color, site, nu, rho);
					}
				}
			}

			//Fifth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[LT::sup(LT::sdn(site,nu),mu)][nu]*this->ridder(unsmearedLattice, site, mu, color, LT::sup(LT::sdn(site,nu),mu), nu, rho);
					}
				}
			}

			//Sixth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[LT::sdn(site,nu)][mu]*this->ridder(unsmearedLattice, site, mu, color, LT::sdn(site,nu), mu, rho);
					}
				}
			}

			//Seventh term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu][color] += smearedDerivative[LT::sdn(site,nu)][nu]*this->ridder(unsmearedLattice, site, mu, color, LT::sdn(site,nu), nu, rho);
					}
				}
			}

		}
	}

	unsmearedDerivative.updateHalo();

	/*//This tensor contains the derivative \frac{\partial U_{ij}}{\partial U_{kl}}
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
	}*/
}



ForceVector SmearingForce::ridder(extended_gauge_lattice_t& unsmearedLattice, int sited, unsigned int mud, int color, int site, unsigned int mu, real_t rho, real_t h) {

	ForceVector a[MAXSIZE][MAXSIZE];
	for (int m = 0; m < MAXSIZE; ++m) {
		for (int n = 0; n < MAXSIZE; ++n) {
			set_to_zero(a[m][n]);
		}
	}

	real_t err = 1000000.;
	real_t factor = REDUCTION*REDUCTION;

	ForceVector result;
	GaugeGroup swap = unsmearedLattice[sited][mud];
	ForceVector tmp = expMap.parameters(unsmearedLattice[sited][mud]);
	tmp[color] += h;
	unsmearedLattice[sited][mud] = expMap.exp(tmp);
	ForceVector plus = expMap.parameters(this->smearLink(unsmearedLattice,site,mu,rho));
	tmp[color] -= 2.*h;
	unsmearedLattice[sited][mud] = expMap.exp(tmp);
	ForceVector minus = expMap.parameters(this->smearLink(unsmearedLattice,site,mu,rho));
	unsmearedLattice[sited][mud] = swap;
	a[0][0] = (plus - minus)/(2.*h);


	result = a[0][0];
	
	for (int m = 1; m < MAXSIZE; ++m) {
		h = h/REDUCTION;
		swap = unsmearedLattice[sited][mud];
		tmp = expMap.parameters(unsmearedLattice[sited][mud]);
		tmp[color] += h;
		unsmearedLattice[sited][mud] = expMap.exp(tmp);
		plus = expMap.parameters(this->smearLink(unsmearedLattice,site,mu,rho));
		tmp[color] -= 2.*h;
		unsmearedLattice[sited][mud] = expMap.exp(tmp);
		minus = expMap.parameters(this->smearLink(unsmearedLattice,site,mu,rho));
		unsmearedLattice[sited][mud] = swap;
		a[0][m] = (plus - minus)/(2.*h);

		factor = REDUCTION*REDUCTION;

		for (int n = 1; n <= m; ++n) {
			a[n][m] = (a[n-1][m]*factor - a[n-1][m-1])/(factor-1.);
			factor = REDUCTION*REDUCTION*factor;
			real_t errt = std::max(fabs(a[n][m][0] - a[n-1][m][0]),fabs(a[n][m][0] - a[n-1][m-1][0]));
			if (errt < err) {
				err = errt;
				result = a[n][m];
			}
		}
	
		if (fabs(a[m][m][0] - a[m-1][m-1][0]) > 2*err) {
			if (err > 0.000000001) {
				std::cout << "Convergence in " << m << " steps and error " << err << " - info "<< sited << " " << mud << "  " << site << " " << mu << std::endl;
			}
			return result;
		}
	}

	if (err > 0.000000001) std::cout << "No convergence in " << MAXSIZE << " steps and error " << err << std::endl;

	
	return result;
}


} /* namespace Update */
