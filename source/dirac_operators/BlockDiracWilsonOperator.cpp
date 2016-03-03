/*
 * BlockDiracWilsonOperator.cpp
 *
 *  Created on: Mar 18, 2013
 *      Author: spiem_01
 */

#include "BlockDiracWilsonOperator.h"
#include "hmc_forces/BlockDiracWilsonFermionForce.h"

namespace Update {

BlockDiracWilsonOperator::BlockDiracWilsonOperator(Color _color) : BlockDiracOperator(_color)  { }

BlockDiracWilsonOperator::BlockDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, Color _color) : BlockDiracOperator(_lattice, _kappa, _color) {
	this->setLattice(_lattice);
}

BlockDiracWilsonOperator::~BlockDiracWilsonOperator() { }

void BlockDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		if (index_lattice[site] == 1) {
			/*GaugeVector tmp[4][4];

			//First we put U(x,mu)*input(x+mu)
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int nu = 0; nu < 4; ++nu) {
					if (index_lattice[DV::sup(site,mu)] == 1) tmp[mu][nu] = lattice[site][mu]*input[DV::sup(site,mu)][nu];
					else set_to_zero(tmp[mu][nu]);
		}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4];
		for (int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] = (tmp[0][0][n]+I*tmp[0][3][n]+tmp[1][0][n]+tmp[1][3][n]+tmp[2][0][n]+I*tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[1][n] = (tmp[0][1][n]+I*tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+tmp[2][1][n]-I*tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
			hopping[2][n] = (I*tmp[0][1][n]-tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+I*tmp[2][0][n]-tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] = (I*tmp[0][0][n]-tmp[0][3][n]-tmp[1][0][n]-tmp[1][3][n]-I*tmp[2][1][n]-tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
		}

		//Then we put U(x-mu,mu)*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				if (index_lattice[LT::sdn(site,mu)] == 1) tmp[mu][nu] = htrans(lattice[LT::sdn(site,mu)][mu])*input[DV::sdn(site,mu)][nu];
				else set_to_zero(tmp[mu][nu]);
		}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] += (tmp[0][0][n]-I*tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+tmp[2][0][n]-I*tmp[2][2][n]+tmp[3][0][n]+tmp[3][2][n]);
			hopping[1][n] += (tmp[0][1][n]-I*tmp[0][2][n]+tmp[1][1][n]+tmp[1][2][n]+tmp[2][1][n]+I*tmp[2][3][n]+tmp[3][1][n]+tmp[3][3][n]);
			hopping[2][n] += (-I*tmp[0][1][n]-tmp[0][2][n]-tmp[1][1][n]-tmp[1][2][n]-I*tmp[2][0][n]-tmp[2][2][n]-tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] += (-I*tmp[0][0][n]-tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+I*tmp[2][1][n]-tmp[2][3][n]-tmp[3][1][n]-tmp[3][3][n]);
		}

		if (gamma5) {
			//The final result is gamma5*input - kappa*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = + input[site][0][n] - kappa*hopping[0][n];
				output[site][1][n] = + input[site][1][n] - kappa*hopping[1][n];
				output[site][2][n] = - (input[site][2][n] + kappa*hopping[2][n]);
				output[site][3][n] = - (input[site][3][n] + kappa*hopping[3][n]);
		}
		}
		else {
			//The final result is input - kappa*gamma5*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = + input[site][0][n] - kappa*hopping[0][n];
				output[site][1][n] = + input[site][1][n] - kappa*hopping[1][n];
				output[site][2][n] = + input[site][2][n] + kappa*hopping[2][n];
				output[site][3][n] = + input[site][3][n] + kappa*hopping[3][n];
		}
		}*/


			//First we start the hopping parameter terms
			GaugeVector tmp_plus[4][2];
			GaugeVector tmp_minus[4][2];

			//Then we project the full spinor in an appropriate half-spinor
			GaugeVector projection_spinor[4][2];

			int site_sup_0 = Vector::sup(site,0);
			int site_sup_1 = Vector::sup(site,1);
			int site_sup_2 = Vector::sup(site,2);
			int site_sup_3 = Vector::sup(site,3);

			for (int n = 0; n < diracVectorLength; ++n) {
				projection_spinor[0][0][n] = std::complex<real_t>(real(input[site_sup_0][0][n])-imag(input[site_sup_0][3][n]),imag(input[site_sup_0][0][n])+real(input[site_sup_0][3][n]));
				projection_spinor[0][1][n] = std::complex<real_t>(real(input[site_sup_0][1][n])-imag(input[site_sup_0][2][n]),imag(input[site_sup_0][1][n])+real(input[site_sup_0][2][n]));
				projection_spinor[1][0][n] = std::complex<real_t>(real(input[site_sup_1][0][n])+real(input[site_sup_1][3][n]),imag(input[site_sup_1][0][n])+imag(input[site_sup_1][3][n]));
				projection_spinor[1][1][n] = std::complex<real_t>(real(input[site_sup_1][1][n])-real(input[site_sup_1][2][n]),imag(input[site_sup_1][1][n])-imag(input[site_sup_1][2][n]));
				projection_spinor[2][0][n] = std::complex<real_t>(real(input[site_sup_2][0][n])-imag(input[site_sup_2][2][n]),imag(input[site_sup_2][0][n])+real(input[site_sup_2][2][n]));
				projection_spinor[2][1][n] = std::complex<real_t>(real(input[site_sup_2][1][n])+imag(input[site_sup_2][3][n]),imag(input[site_sup_2][1][n])-real(input[site_sup_2][3][n]));
				projection_spinor[3][0][n] = std::complex<real_t>(real(input[site_sup_3][0][n])-real(input[site_sup_3][2][n]),imag(input[site_sup_3][0][n])-imag(input[site_sup_3][2][n]));
				projection_spinor[3][1][n] = std::complex<real_t>(real(input[site_sup_3][1][n])-real(input[site_sup_3][3][n]),imag(input[site_sup_3][1][n])-imag(input[site_sup_3][3][n]));
			}

			//Now we can put U(x,mu)*input(x+mu)
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int nu = 0; nu < 2; ++nu) {
					if (index_lattice[Lattice::sup(site,mu)] == 1) tmp_plus[mu][nu] = lattice[site][mu]*projection_spinor[mu][nu];
					else set_to_zero(tmp_plus[mu][nu]);
				}
			}

			int site_down_0 = Vector::sdn(site,0);
			int site_down_1 = Vector::sdn(site,1);
			int site_down_2 = Vector::sdn(site,2);
			int site_down_3 = Vector::sdn(site,3);

			for (int n = 0; n < diracVectorLength; ++n) {
				projection_spinor[0][0][n] = std::complex<real_t>(real(input[site_down_0][0][n])+imag(input[site_down_0][3][n]),imag(input[site_down_0][0][n])-real(input[site_down_0][3][n]));
				projection_spinor[0][1][n] = std::complex<real_t>(real(input[site_down_0][1][n])+imag(input[site_down_0][2][n]),imag(input[site_down_0][1][n])-real(input[site_down_0][2][n]));
				projection_spinor[1][0][n] = std::complex<real_t>(real(input[site_down_1][0][n])-real(input[site_down_1][3][n]),imag(input[site_down_1][0][n])-imag(input[site_down_1][3][n]));
				projection_spinor[1][1][n] = std::complex<real_t>(real(input[site_down_1][1][n])+real(input[site_down_1][2][n]),imag(input[site_down_1][1][n])+imag(input[site_down_1][2][n]));
				projection_spinor[2][0][n] = std::complex<real_t>(real(input[site_down_2][0][n])+imag(input[site_down_2][2][n]),imag(input[site_down_2][0][n])-real(input[site_down_2][2][n]));
				projection_spinor[2][1][n] = std::complex<real_t>(real(input[site_down_2][1][n])-imag(input[site_down_2][3][n]),imag(input[site_down_2][1][n])+real(input[site_down_2][3][n]));
				projection_spinor[3][0][n] = std::complex<real_t>(real(input[site_down_3][0][n])+real(input[site_down_3][2][n]),imag(input[site_down_3][0][n])+imag(input[site_down_3][2][n]));
				projection_spinor[3][1][n] = std::complex<real_t>(real(input[site_down_3][1][n])+real(input[site_down_3][3][n]),imag(input[site_down_3][1][n])+imag(input[site_down_3][3][n]));
			}

			//Then we put U(x-mu,mu)*input(x-mu)
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int nu = 0; nu < 2; ++nu) {
					if (index_lattice[Lattice::sdn(site,mu)] == 1) {
						GaugeVector tmp = htrans(lattice[Lattice::sdn(site,mu)][mu])*projection_spinor[mu][nu];
						tmp_minus[mu][nu] = tmp_plus[mu][nu] - tmp;
						tmp_plus[mu][nu] += tmp;
					}
					else {
						tmp_minus[mu][nu] = tmp_plus[mu][nu];
					}
				}
			}


			if (gamma5) {
				//The final result is gamma5*input - kappa*hopping
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = std::complex<real_t>(0.,twist)*input[site][0][n] + std::complex<real_t>(real(input[site][0][n])-kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n]) -kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
					output[site][1][n] = std::complex<real_t>(0.,twist)*input[site][1][n] + std::complex<real_t>(real(input[site][1][n])-kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n]) -kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
					output[site][2][n] = std::complex<real_t>(0.,twist)*input[site][2][n] + std::complex<real_t>(-(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n]))),-(imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n]))));
					output[site][3][n] = std::complex<real_t>(0.,twist)*input[site][3][n] + std::complex<real_t>(-(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n]))),-(imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n]))));
				}
			}
			else {
				//The final result is input - kappa*gamma5*hopping
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = std::complex<real_t>(0.,twist)*input[site][0][n] + std::complex<real_t>(real(input[site][0][n]) - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n]) - kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
					output[site][1][n] = std::complex<real_t>(0.,twist)*input[site][1][n] + std::complex<real_t>(real(input[site][1][n]) - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n]) - kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
					output[site][2][n] = std::complex<real_t>(0.,-twist)*input[site][2][n] + std::complex<real_t>(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])),imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
					output[site][3][n] = std::complex<real_t>(0.,-twist)*input[site][3][n] + std::complex<real_t>(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])),imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
				}
			}
		}
		else {
			if (gamma5) {
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = std::complex<real_t>(1.,twist)*input[site][0][n];
					output[site][1][n] = std::complex<real_t>(1.,twist)*input[site][1][n];
					output[site][2][n] = std::complex<real_t>(-1.,twist)*input[site][2][n];
					output[site][3][n] = std::complex<real_t>(-1.,twist)*input[site][3][n];
				}
			}
			else {
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = + std::complex<real_t>(1.,twist)*input[site][0][n];
					output[site][1][n] = + std::complex<real_t>(1.,twist)*input[site][1][n];
					output[site][2][n] = + std::complex<real_t>(1.,-twist)*input[site][2][n];
					output[site][3][n] = + std::complex<real_t>(1.,-twist)*input[site][3][n];
				}
			}
		}
	}
	//output.updateHalo(); Update halo called in the outer steps!
}

void BlockDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const std::complex<real_t>& alpha) {
	typedef reduced_fermion_lattice_t LT;
	typedef reduced_dirac_vector_t DV;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		if (index_lattice[site] == 1) {
			//First we start the hopping parameter terms
			GaugeVector tmp[4][4];

			//First we put U(x,mu)*input(x+mu)
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int nu = 0; nu < 4; ++nu) {
					tmp[mu][nu] = lattice[site][mu]*vector1[DV::sup(site,mu)][nu];
				}
			}

			//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
			GaugeVector hopping[4];
			for (int n = 0; n < diracVectorLength; ++n) {
				hopping[0][n] = (tmp[0][0][n]+I*tmp[0][3][n]+tmp[1][0][n]+tmp[1][3][n]+tmp[2][0][n]+I*tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
				hopping[1][n] = (tmp[0][1][n]+I*tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+tmp[2][1][n]-I*tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
				hopping[2][n] = (I*tmp[0][1][n]-tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+I*tmp[2][0][n]-tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
				hopping[3][n] = (I*tmp[0][0][n]-tmp[0][3][n]-tmp[1][0][n]-tmp[1][3][n]-I*tmp[2][1][n]-tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
			}

			//Then we put U(x-mu,mu)*input(x-mu)
			for (unsigned int mu = 0; mu < 4; ++mu) {
				for (unsigned int nu = 0; nu < 4; ++nu) {
					tmp[mu][nu] = htrans(lattice[LT::sdn(site,mu)][mu])*vector1[DV::sdn(site,mu)][nu];
				}
			}

			//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
			for (int n = 0; n < diracVectorLength; ++n) {
				hopping[0][n] += (tmp[0][0][n]-I*tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+tmp[2][0][n]-I*tmp[2][2][n]+tmp[3][0][n]+tmp[3][2][n]);
				hopping[1][n] += (tmp[0][1][n]-I*tmp[0][2][n]+tmp[1][1][n]+tmp[1][2][n]+tmp[2][1][n]+I*tmp[2][3][n]+tmp[3][1][n]+tmp[3][3][n]);
				hopping[2][n] += (-I*tmp[0][1][n]-tmp[0][2][n]-tmp[1][1][n]-tmp[1][2][n]-I*tmp[2][0][n]-tmp[2][2][n]-tmp[3][0][n]-tmp[3][2][n]);
				hopping[3][n] += (-I*tmp[0][0][n]-tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+I*tmp[2][1][n]-tmp[2][3][n]-tmp[3][1][n]-tmp[3][3][n]);
			}

			if (gamma5) {
				//The final result is gamma5*input - kappa*hopping
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(1.,twist)*vector1[site][0][n] - kappa*hopping[0][n];
					output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(1.,twist)*vector1[site][1][n] - kappa*hopping[1][n];
					output[site][2][n] = alpha*vector2[site][2][n] - (std::complex<real_t>(1.,-twist)*vector1[site][2][n] + kappa*hopping[2][n]);
					output[site][3][n] = alpha*vector2[site][3][n] - (std::complex<real_t>(1.,-twist)*vector1[site][3][n] + kappa*hopping[3][n]);
				}
			}
			else {
				//The final result is input - kappa*gamma5*hopping
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(1.,twist)*vector1[site][0][n] - kappa*hopping[0][n];
					output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(1.,twist)*vector1[site][1][n] - kappa*hopping[1][n];
					output[site][2][n] = alpha*vector2[site][2][n] + std::complex<real_t>(1.,-twist)*vector1[site][2][n] + kappa*hopping[2][n];
					output[site][3][n] = alpha*vector2[site][3][n] + std::complex<real_t>(1.,-twist)*vector1[site][3][n] + kappa*hopping[3][n];
				}
			}
		}
		else {
			if (gamma5) {
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(1.,twist)*vector1[site][0][n];
					output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(1.,twist)*vector1[site][1][n];
					output[site][2][n] = alpha*vector2[site][2][n] + std::complex<real_t>(-1.,twist)*vector1[site][2][n];
					output[site][3][n] = alpha*vector2[site][3][n] + std::complex<real_t>(-1.,twist)*vector1[site][3][n];
				}
			}
			else {
				for (int n = 0; n < diracVectorLength; ++n) {
					output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(1.,twist)*vector1[site][0][n];
					output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(1.,twist)*vector1[site][1][n];
					output[site][2][n] = alpha*vector2[site][2][n] + std::complex<real_t>(1.,-twist)*vector1[site][2][n];
					output[site][3][n] = alpha*vector2[site][3][n] + std::complex<real_t>(1.,-twist)*vector1[site][3][n];
				}
			}
		}
	}
	//output.updateHalo();//TODO is needed?
}

FermionForce* BlockDiracWilsonOperator::getForce() const {
	return new BlockDiracWilsonFermionForce(kappa, xBlockSize, yBlockSize, zBlockSize, tBlockSize);
}

} /* namespace Update */
