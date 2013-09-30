/*
 * ImprovedFermionForce.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: spiem_01
 */

#include "ImprovedFermionForce.h"
#include "Gamma.h"

namespace Update {

ImprovedFermionForce::ImprovedFermionForce(real_t _kappa, real_t _csw) : DiracWilsonFermionForce(_kappa), csw(_csw) { }

ImprovedFermionForce::~ImprovedFermionForce() { }

FermionicForceMatrix ImprovedFermionForce::derivative(const extended_fermion_lattice_t& lattice, const extended_dirac_vector_t& X, const extended_dirac_vector_t& Y, int site, int mu) const {
	typedef extended_dirac_vector_t Vector;
	typedef extended_fermion_lattice_t Lattice;
	FermionicForceMatrix force;
	set_to_zero(force);
	//First Term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[site][beta];
				}
				force += tensor(X[site][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[site][beta];
				}
				force += tensor(Y[site][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu]);
			FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu])*htrans(lattice[site][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[site][beta];
				}
				force += tensor(linkMatrixRigth*X[site][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[site][beta];
				}
				force += tensor(linkMatrixRigth*Y[site][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Second term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu]);
			FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu])*lattice[Lattice::sdn(site,nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[site][beta];
				}
				force -= tensor(linkMatrixRigth*X[site][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[site][beta];
				}
				force -= tensor(linkMatrixRigth*Y[site][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu])*lattice[Lattice::sdn(site,nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[site][beta];
				}
				force -= tensor(X[site][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[site][beta];
				}
				force -= tensor(Y[site][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Third term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,nu)][beta];
				}
				force += tensor(lattice[site][nu]*X[Vector::sup(site,nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,nu)][beta];
				}
				force += tensor(lattice[site][nu]*Y[Vector::sup(site,nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu];
			FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu]*htrans(lattice[Lattice::sup(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,nu)][beta];
				}
				force += tensor(linkMatrixRigth*X[Vector::sup(site,nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,nu)][beta];
				}
				force += tensor(linkMatrixRigth*Y[Vector::sup(site,nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Fourth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu]);
			FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(site,nu)][beta];
				}
				force -= tensor(linkMatrixRigth*X[Vector::sdn(site,nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(site,nu)][beta];
				}
				force -= tensor(linkMatrixRigth*Y[Vector::sdn(site,nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu])*htrans(lattice[Lattice::sdn(site,nu)][mu]);
			FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(site,nu)][beta];
				}
				force -= tensor(linkMatrixRigth*X[Vector::sdn(site,nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(site,nu)][beta];
				}
				force -= tensor(linkMatrixRigth*Y[Vector::sdn(site,nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Fifth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicGroup linkMatrixRigth = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu]*htrans(lattice[Lattice::sup(site,mu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
				}
				force += tensor(linkMatrixRigth*X[Vector::sup(site,mu)][alpha],(I*kappa*csw)*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,mu)][beta];
				}
				force += tensor(linkMatrixRigth*Y[Vector::sup(site,mu)][alpha],(I*kappa*csw)*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu]*htrans(lattice[Lattice::sup(site,mu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
				}
				force += tensor(lattice[site][mu]*X[Vector::sup(site,mu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,mu)][beta];
				}
				force += tensor(lattice[site][mu]*Y[Vector::sup(site,mu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Sixth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu]*lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
				}
				force -= tensor(lattice[site][mu]*X[Vector::sup(site,mu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,mu)][beta];
				}
				force -= tensor(lattice[site][mu]*Y[Vector::sup(site,mu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu]*lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(site,mu)][beta];
				}
				force -= tensor(linkMatrixRigth*X[Vector::sup(site,mu)][alpha],(I*kappa*csw)*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(site,mu)][beta];
				}
				force -= tensor(linkMatrixRigth*Y[Vector::sup(site,mu)][alpha],(I*kappa*csw)*tmp);
			}
		}
	}

	//Seventh term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*lattice[Lattice::sup(site,mu)][nu];
			FermionicGroup linkMatrixRigth = lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(Vector::sup(site,mu),nu)][beta];
				}
				force += tensor(linkMatrixRigth*X[Vector::sup(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(Vector::sup(site,mu),nu)][beta];
				}
				force += tensor(linkMatrixRigth*Y[Vector::sup(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*lattice[site][nu]*lattice[Lattice::sup(site,nu)][mu];
			FermionicGroup linkMatrixRigth = lattice[site][mu]*lattice[Lattice::sup(site,mu)][nu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sup(Vector::sup(site,mu),nu)][beta];
				}
				force += tensor(linkMatrixRigth*X[Vector::sup(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sup(Vector::sup(site,mu),nu)][beta];
				}
				force += tensor(linkMatrixRigth*Y[Vector::sup(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	//Eighth term
	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[site][mu])*htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu];
			FermionicGroup linkMatrixRigth = lattice[site][mu]*htrans(lattice[Lattice::sdn(Vector::sup(site,mu),nu)][nu]);
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(Vector::sup(site,mu),nu)][beta];
				}
				force -= tensor(linkMatrixRigth*X[Vector::sdn(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(Vector::sup(site,mu),nu)][beta];
				}
				force -= tensor(linkMatrixRigth*Y[Vector::sdn(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	for (int nu = 0; nu < 4; ++nu) {
		if (nu != mu) {
			FermionicForceMatrix linkMatrixLeft = (I*kappa*csw)*htrans(lattice[Lattice::sdn(Lattice::sup(site,mu),nu)][nu]);
			FermionicGroup linkMatrixRigth = htrans(lattice[Lattice::sdn(site,nu)][nu])*lattice[Lattice::sdn(site,nu)][mu];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				GaugeVector tmp;
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*Y[Vector::sdn(Vector::sup(site,mu),nu)][beta];
				}
				force -= tensor(linkMatrixRigth*X[Vector::sdn(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
				//Now we exchange x with y
				set_to_zero(tmp);
				for (unsigned int beta = 0; beta < 4; ++beta) {
					if (Sigma::g5sigma(mu,nu,alpha,beta) != 0.) tmp += Sigma::g5sigma(mu,nu,alpha,beta)*X[Vector::sdn(Vector::sup(site,mu),nu)][beta];
				}
				force -= tensor(linkMatrixRigth*Y[Vector::sdn(Vector::sup(site,mu),nu)][alpha],linkMatrixLeft*tmp);
			}
		}
	}

	return (force/8.) + DiracWilsonFermionForce::derivative(lattice, X, Y, site, mu);
}

} /* namespace Update */
