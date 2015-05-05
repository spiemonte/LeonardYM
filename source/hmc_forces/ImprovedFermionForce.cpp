/*
 * ImprovedFermionForce.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: spiem_01
 */

#include "ImprovedFermionForce.h"
#include "utils/Gamma.h"

namespace Update {

ImprovedFermionForce::ImprovedFermionForce(real_t _kappa, real_t _csw) : DiracWilsonFermionForce(_kappa), csw(_csw) { }

ImprovedFermionForce::ImprovedFermionForce(const ImprovedFermionForce& toCopy) : DiracWilsonFermionForce(toCopy.kappa), csw(toCopy.csw) { }

ImprovedFermionForce::~ImprovedFermionForce() { }

FermionicForceMatrix ImprovedFermionForce::derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const {
	typedef extended_dirac_vector_t Vector;
	typedef extended_fermion_lattice_t Lattice;
	FermionicForceMatrix force;
	set_to_zero(force);

	std::complex<real_t> ikc = (I*kappa*csw);

	GaugeVector projSpinorX[4][4];
	GaugeVector projSpinorY[4][4];

	//We project the spinors for the first two terms
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[site][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[site][beta];
					}
				}
			}
		}
	}

	//First Term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(X[site][alpha],latticeForcesLeft[0](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(Y[site][alpha],latticeForcesLeft[0](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu]);
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[1](site,mu,nu)*X[site][alpha],latticeForcesLeft[1](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[1](site,mu,nu)*Y[site][alpha],latticeForcesLeft[1](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/

	//Second term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu]);
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu])*lattice[Lattice::sdn(site,nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[1](site,mu,nu)*X[site][alpha],latticeForcesLeft[1](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[1](site,mu,nu)*Y[site][alpha],latticeForcesLeft[1](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu])*lattice[Lattice::sdn(site,nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(X[site][alpha],latticeForcesLeft[3](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(Y[site][alpha],latticeForcesLeft[3](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/

	//We project the spinors for the next term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,nu)][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,nu)][beta];
					}
				}
			}
		}
	}

	//Third term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[2](site,mu,nu)*X[Vector::sup(site,nu)][alpha],latticeForcesLeft[2](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[2](site,mu,nu)*Y[Vector::sup(site,nu)][alpha],latticeForcesLeft[2](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu];
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[5](site,mu,nu)*X[Vector::sup(site,nu)][alpha],latticeForcesLeft[5](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[5](site,mu,nu)*Y[Vector::sup(site,nu)][alpha],latticeForcesLeft[5](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/


	//We project the spinors for the next term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(site,nu)][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(site,nu)][beta];
					}
				}
			}
		}
	}

	//Fourth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu]);
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[3](site,mu,nu)*X[Vector::sdn(site,nu)][alpha],latticeForcesLeft[3](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[3](site,mu,nu)*Y[Vector::sdn(site,nu)][alpha],latticeForcesLeft[3](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu]);
			//FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[7](site,mu,nu)*X[Vector::sdn(site,nu)][alpha],latticeForcesLeft[7](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[7](site,mu,nu)*Y[Vector::sdn(site,nu)][alpha],latticeForcesLeft[7](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/


	//We project the spinors for the next two terms
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,mu)][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
					}
				}
			}
		}
	}

	//Fifth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicGroup linkMatrixRigth = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu]*htrans(lattice[Lattice::sup(site,mu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[4](site,mu,nu)*X[Vector::sup(site,mu)][alpha],projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[4](site,mu,nu)*Y[Vector::sup(site,mu)][alpha],projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu]*htrans(lattice[Lattice::sup(site,mu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[9](site,mu,nu)*X[Vector::sup(site,mu)][alpha],latticeForcesLeft[9](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[9](site,mu,nu)*Y[Vector::sup(site,mu)][alpha],latticeForcesLeft[9](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/

	//Sixth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu]*lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[5](site,mu,nu)*X[Vector::sup(site,mu)][alpha],latticeForcesLeft[5](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[5](site,mu,nu)*Y[Vector::sup(site,mu)][alpha],latticeForcesLeft[5](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu]*lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[11](site,mu,nu)*X[Vector::sup(site,mu)][alpha],projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[11](site,mu,nu)*Y[Vector::sup(site,mu)][alpha],projSpinorX[nu][alpha]);
			}
		}
	}*/


	//We project the spinors for the next term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(Vector::sup(site,mu),nu)][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(Vector::sup(site,mu),nu)][beta];
					}
				}
			}
		}
	}

	//Seventh term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu];
			//FermionicGroup linkMatrixRigth = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[6](site,mu,nu)*X[Vector::sup(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[6](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[6](site,mu,nu)*Y[Vector::sup(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[6](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu];
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force += tensor(latticeForcesRight[13](site,mu,nu)*X[Vector::sup(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[13](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force += tensor(latticeForcesRight[13](site,mu,nu)*Y[Vector::sup(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[13](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/

	//We project the spinors for the next term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				set_to_zero(projSpinorX[nu][alpha]);
				set_to_zero(projSpinorY[nu][alpha]);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != static_cast<real_t>(0.)) {
						projSpinorX[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(Vector::sup(site,mu),nu)][beta];
						projSpinorY[nu][alpha] += ikc*Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(Vector::sup(site,mu),nu)][beta];
					}
				}
			}
		}
	}

	//Eighth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu];
			//FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Vector::sup(site,mu),nu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[7](site,mu,nu)*X[Vector::sdn(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[7](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[7](site,mu,nu)*Y[Vector::sdn(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[7](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}

	/*for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			//FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu]);
			//FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				force -= tensor(latticeForcesRight[15](site,mu,nu)*X[Vector::sdn(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[15](site,mu,nu)*projSpinorY[nu][alpha]);
				//Now we exchange x with y
				force -= tensor(latticeForcesRight[15](site,mu,nu)*Y[Vector::sdn(Vector::sup(site,mu),nu)][alpha],latticeForcesLeft[15](site,mu,nu)*projSpinorX[nu][alpha]);
			}
		}
	}*/

	return (force/4.) + DiracWilsonFermionForce::derivative(lattice, X, Y, site, mu);
}

void ImprovedFermionForce::setLattice(const extended_fermion_lattice_t& lattice) {
	typedef extended_fermion_lattice_t Lattice;

	//We allocate the memory
	for (int i = 0; i < 8; ++i) {
		latticeForcesLeft[i].setLocaLatticeSize(lattice.localsize);
		latticeForcesRight[i].setLocaLatticeSize(lattice.localsize);
	}

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		for (int mu = 0; mu < 4; ++mu) {
			for (int nu = 0; nu < 4; ++nu) {
				if (nu != mu) {
					//First Term
					latticeForcesLeft[0](site,mu,nu) = lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu])*htrans(lattice[site][nu]);

					//Second term
					latticeForcesLeft[1](site,mu,nu) = htrans(lattice[site][mu]);
					latticeForcesRight[1](site,mu,nu) = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu])*lattice[Lattice::sdn(site,nu)][nu];

					//Third Term
					latticeForcesLeft[2](site,mu,nu) = lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu]);
					latticeForcesRight[2](site,mu,nu) = lattice[site][nu];

					//Fourth term
					latticeForcesLeft[3](site,mu,nu) = htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu]);
					latticeForcesRight[3](site,mu,nu) = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu]);

					//Fifth term
					latticeForcesRight[4](site,mu,nu) = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu]*htrans(lattice[Lattice::sup(site,mu)][nu]);

					//Sixth term
					latticeForcesLeft[5](site,mu,nu) = htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu]*lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu];
					latticeForcesRight[5](site,mu,nu) = lattice[site][mu];

					//Seventh term
					latticeForcesLeft[6](site,mu,nu) = lattice[Lattice::sup(site,mu)][nu];
					latticeForcesRight[6](site,mu,nu) = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu];

					//Eighth term
					latticeForcesLeft[7](site,mu,nu) = htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu];
					latticeForcesRight[7](site,mu,nu) = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu]);
				}
			}
		}
	}
}

} /* namespace Update */
