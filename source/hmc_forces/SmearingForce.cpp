#include "SmearingForce.h"
#define REDUCTION 1.3
#define MAXSIZE 25

namespace Update {

SmearingForce::SmearingForce() : StoutSmearing(), expMap() { }

SmearingForce::~SmearingForce() { }

void SmearingForce::force(const extended_fermion_force_lattice_t& actionDerivative, const extended_gauge_lattice_t& originalUnsmearedLattice, extended_gauge_lattice_t& unsmearedDerivative, real_t rho) {
	typedef extended_gauge_lattice_t LT;
	
#ifndef MULTITHREADING
	extended_gauge_lattice_t unsmearedLattice = originalUnsmearedLattice;
	for (int site = 0; site < unsmearedLattice.localsize; ++site) {
#endif
#ifdef MULTITHREADING
	extended_gauge_lattice_t thUnsmearedLattice[omp_get_max_threads()];
	for (int i = 0; i < omp_get_max_threads(); ++i) thUnsmearedLattice[i] = originalUnsmearedLattice;
#pragma omp parallel for
	for (int site = 0; site < originalUnsmearedLattice.localsize; ++site) {
		extended_gauge_lattice_t& unsmearedLattice = thUnsmearedLattice[omp_get_thread_num()];
#endif
		for (unsigned int mu = 0; mu < 4; ++mu) {
			set_to_zero(unsmearedDerivative[site][mu]);
			//First term
			for (int color = 0; color < numberColors*numberColors - 1; ++color) {
				unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, site, mu, rho))*gaugeLieGenerators.get(color);
			}

			//Second term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, LT::sup(site,mu), nu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

			//Third term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, LT::sup(site,nu), mu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

			//Fourth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, site, nu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

			//Fifth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, LT::sup(LT::sdn(site,nu),mu), nu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

			//Sixth term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, LT::sdn(site,nu), mu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

			//Seventh term
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					for (int color = 0; color < numberColors*numberColors - 1; ++color) {
						unsmearedDerivative[site][mu] += std::complex<real_t>(0,this->ridder(actionDerivative, unsmearedLattice, site, mu, color, LT::sdn(site,nu), nu, rho))*gaugeLieGenerators.get(color);
					}
				}
			}

		}
	}

	unsmearedDerivative.updateHalo();
}



real_t SmearingForce::ridder(const extended_fermion_force_lattice_t& actionDerivative, extended_gauge_lattice_t& unsmearedLattice, int sited, unsigned int mud, int color, int site, unsigned int mu, real_t rho, real_t h) {
	real_t a[MAXSIZE][MAXSIZE];
	for (int m = 0; m < MAXSIZE; ++m) {
		for (int n = 0; n < MAXSIZE; ++n) {
			a[m][n] = 0;
		}
	}

	real_t err = 1000000.;
	real_t factor = REDUCTION*REDUCTION;

	GaugeGroup swap = unsmearedLattice[sited][mud];
	
	ForceVector tmp;
	set_to_zero(tmp);
	tmp[color] = h;
	expMap.exp(tmp, unsmearedLattice[sited][mud]);
	unsmearedLattice[sited][mud] = unsmearedLattice[sited][mud]*swap;
#ifdef ADJOINT
	FermionicGroup link;
	ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(this->smearLink(unsmearedLattice,site,mu,rho),link);
	real_t plus = real(trace(actionDerivative[site][mu]*link));
#endif
#ifndef ADJOINT
	real_t plus = real(trace(actionDerivative[site][mu]*this->smearLink(unsmearedLattice,site,mu,rho)));
#endif
	tmp[color] = -h;
	expMap.exp(tmp, unsmearedLattice[sited][mud]);
	unsmearedLattice[sited][mud] = unsmearedLattice[sited][mud]*swap;
#ifdef ADJOINT
	ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(this->smearLink(unsmearedLattice,site,mu,rho),link);
	real_t minus = real(trace(actionDerivative[site][mu]*link));
#endif
#ifndef ADJOINT
	real_t minus = real(trace(actionDerivative[site][mu]*this->smearLink(unsmearedLattice,site,mu,rho)));
#endif
	unsmearedLattice[sited][mud] = swap;
	
	a[0][0] = (plus - minus)/(2.*h);


	real_t result = a[0][0];
	
	for (int m = 1; m < MAXSIZE; ++m) {
		h = h/REDUCTION;
		
		swap = unsmearedLattice[sited][mud];
		
		tmp[color] = h;
		expMap.exp(tmp, unsmearedLattice[sited][mud]);
		unsmearedLattice[sited][mud] = unsmearedLattice[sited][mud]*swap;
#ifdef ADJOINT
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(this->smearLink(unsmearedLattice,site,mu,rho),link);
		plus = real(trace(actionDerivative[site][mu]*link));
#endif
#ifndef ADJOINT
		plus = real(trace(actionDerivative[site][mu]*this->smearLink(unsmearedLattice,site,mu,rho)));
#endif
		tmp[color] = -h;
		expMap.exp(tmp, unsmearedLattice[sited][mud]);
		unsmearedLattice[sited][mud] = unsmearedLattice[sited][mud]*swap;
#ifdef ADJOINT
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(this->smearLink(unsmearedLattice,site,mu,rho),link);
		minus = real(trace(actionDerivative[site][mu]*link));
#endif
#ifndef ADJOINT
		minus = real(trace(actionDerivative[site][mu]*this->smearLink(unsmearedLattice,site,mu,rho)));
#endif
		unsmearedLattice[sited][mud] = swap;
		a[0][m] = (plus - minus)/(2.*h);

		factor = REDUCTION*REDUCTION;

		for (int n = 1; n <= m; ++n) {
			a[n][m] = (a[n-1][m]*factor - a[n-1][m-1])/(factor-1.);
			factor = REDUCTION*REDUCTION*factor;
			real_t errt = std::max(std::abs(a[n][m] - a[n-1][m]),std::abs(a[n][m] - a[n-1][m-1]));
			if (errt < err) {
				err = errt;
				result = a[n][m];
			}
		}
	
		if (std::abs(a[m][m] - a[m-1][m-1]) > 2*err) {
			if (err > 0.0000001) {
				std::cout << "SmearingForce::Bad convergence in " << m << " steps and error " << err << std::endl;
			}
			return result;
		}
	}

	if (err > 0.0000001) std::cout << "SmearingForce::Bad convergence in " << MAXSIZE << " steps and error " << err << std::endl;

	
	return result;
}


} /* namespace Update */
