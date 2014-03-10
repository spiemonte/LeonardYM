/*
 * ImprovedDiracWilsonOperator.cpp
 *
 *  Created on: May 4, 2012
 *      Author: spiem_01
 */

#include "ImprovedDiracWilsonOperator.h"
#include "../ImprovedFermionForce.h"

namespace Update {

inline real_t conj(const real_t& t) {
	return t;
}

inline std::complex<real_t> multiply_by_I(const std::complex<real_t>& a) {
	return std::complex<real_t>(-imag(a),real(a));
}

ImprovedDiracWilsonOperator::ImprovedDiracWilsonOperator() : DiracOperator(), csw(0.) { }

ImprovedDiracWilsonOperator::ImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, double _kappa, double _csw, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5), csw(_csw) {
	this->updateFieldStrength(_lattice);
}

ImprovedDiracWilsonOperator::~ImprovedDiracWilsonOperator() { }

void ImprovedDiracWilsonOperator::multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;
	const reduced_fermion_lattice_t& linkconf = (lattice);

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		//the best
		std::complex<real_t> projection_spinor_minus[diracVectorLength], projection_spinor_plus[diracVectorLength], tmm, tmp;
		{
			{
				const size_t site_down = Vector::sdn(site,0);
				const size_t site_up = Vector::sup(site,0);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(input[site_down][0][n].real()+input[site_down][3][n].imag(),input[site_down][0][n].imag()-input[site_down][3][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(input[site_up][0][n].real()-input[site_up][3][n].imag(),input[site_up][0][n].imag()+input[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][0](i,n);
					}
					output[site][0][i] = input[site][0][i] - kappa*(tmp+tmm);
					output[site][3][i] = -input[site][3][i] + kappa*std::complex<real_t>(tmm.imag() - tmp.imag(),tmp.real() - tmm.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(input[site_down][1][n].real()+input[site_down][2][n].imag(), input[site_down][1][n].imag()-input[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(input[site_up][1][n].real()-input[site_up][2][n].imag(),input[site_up][1][n].imag()+input[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][0](i,n);
					}

					output[site][1][i] = input[site][1][i]- kappa*(tmp+tmm);
					output[site][2][i] = -input[site][2][i]+ kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real()-tmm.real());//(positI(tmp) +negatI(tmm));
				}
			}
			{
				const size_t site_down = Vector::sdn(site,1);
				const size_t site_up = Vector::sup(site,1);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = input[site_down][0][n] - (input[site_down][3][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = input[site_up][0][n] + (input[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][1](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*(tmm-tmp);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = input[site_down][1][n] + (input[site_down][2][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = input[site_up][1][n] - (input[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][1](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*(tmp-tmm);
				}
			}
			{
				const size_t site_down = Vector::sdn(site,2);
				const size_t site_up = Vector::sup(site,2);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(input[site_down][0][n].real() + input[site_down][2][n].imag(), input[site_down][0][n].imag() - input[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(input[site_up][0][n].real() - input[site_up][2][n].imag(), input[site_up][0][n].imag() + input[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][2](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real() - tmm.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(input[site_down][1][n].real() - input[site_down][3][n].imag(), input[site_down][1][n].imag() + input[site_down][3][n].real());
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(input[site_up][1][n].real() + input[site_up][3][n].imag(), input[site_up][1][n].imag() - input[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][2](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*std::complex<real_t>(tmp.imag() - tmm.imag(), tmm.real() - tmp.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,3);
				const size_t site_up = Vector::sup(site,3);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = input[site_down][0][n] + (input[site_down][2][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = input[site_up][0][n] - (input[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][3](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*(tmp-tmm);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = input[site_down][1][n] + (input[site_down][3][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = input[site_up][1][n] - (input[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][3](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*(tmp-tmm);
				}
			}

			//We store the result of the clover term in an intermediate vector
			GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				clover[0][i] = 0;
				clover[1][i] = 0;
				clover[2][i] = 0;
				clover[3][i] = 0;
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += multiply_by_I((-F[site][0].at(i,j)+F[site][5].at(i,j))*input[site][0][j]);
					clover[1][i] += multiply_by_I((+F[site][0].at(i,j)-F[site][5].at(i,j))*input[site][1][j]);
					clover[2][i] += multiply_by_I((+F[site][0].at(i,j)+F[site][5].at(i,j))*input[site][2][j]);
					clover[3][i] += multiply_by_I((-F[site][0].at(i,j)-F[site][5].at(i,j))*input[site][3][j]);
				}
			}

#ifdef ADJOINT
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += std::complex<real_t>(+F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*input[site][1][j];
					clover[1][i] += std::complex<real_t>(-F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*input[site][0][j];
					clover[2][i] += std::complex<real_t>(-F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*input[site][3][j];
					clover[3][i] += std::complex<real_t>(+F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*input[site][2][j];
				}
			}
#endif
#ifndef ADJOINT
			//We store the result of the clover term in an intermediate vector
			//GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += (+F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*input[site][1][j];
					clover[1][i] += (-F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*input[site][0][j];
					clover[2][i] += (-F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*input[site][3][j];
					clover[3][i] += (+F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*input[site][2][j];
				}
			}
#endif

			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] += (kappa*csw)*clover[0][n];
				output[site][1][n] += (kappa*csw)*clover[1][n];
				output[site][2][n] += (kappa*csw)*clover[2][n];
				output[site][3][n] += (kappa*csw)*clover[3][n];
			}
		}
		if (!gamma5) {
			for (int i = 0; i < diracVectorLength; ++i) {
				output[site][2][i] = -output[site][2][i];
				output[site][3][i] = -output[site][3][i];
			}
		}

/*
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
				tmp_plus[mu][nu] = lattice[site][mu]*projection_spinor[mu][nu];
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
				GaugeVector tmp = htrans(lattice[Lattice::sdn(site,mu)][mu])*projection_spinor[mu][nu];
				tmp_minus[mu][nu] = tmp_plus[mu][nu] - tmp;
				tmp_plus[mu][nu] += tmp;
			}
		}

		//Clover matrix needed
		std::complex<real_t> FA[4][diracVectorLength][diracVectorLength];//TODO why real?
		for (int i = 0; i < diracVectorLength; ++i) {
			for (int j = 0; j < diracVectorLength; ++j) {
				FA[0][i][j] = -F[site][0].at(i,j)+F[site][5].at(i,j);
				FA[1][i][j] = +F[site][0].at(i,j)-F[site][5].at(i,j);
				FA[2][i][j] = +F[site][0].at(i,j)+F[site][5].at(i,j);
				FA[3][i][j] = -F[site][0].at(i,j)-F[site][5].at(i,j);
			}
		}

		std::complex<real_t> FB[4][diracVectorLength][diracVectorLength];
		for (int i = 0; i < diracVectorLength; ++i) {
			for (int j = 0; j < diracVectorLength; ++j) {
				FB[0][i][j] = +F[site][1].at(i,j)+F[site][4].at(i,j)+I*(+F[site][2].at(i,j)-F[site][3].at(i,j));
				FB[1][i][j] = -F[site][1].at(i,j)-F[site][4].at(i,j)+I*(+F[site][2].at(i,j)-F[site][3].at(i,j));
				FB[2][i][j] = -F[site][1].at(i,j)+F[site][4].at(i,j)+I*(+F[site][2].at(i,j)+F[site][3].at(i,j));
				FB[3][i][j] = +F[site][1].at(i,j)-F[site][4].at(i,j)+I*(+F[site][2].at(i,j)+F[site][3].at(i,j));
			}
		}

		//We store the result of the clover term in an intermediate vector
		GaugeVector clover[4];
		for (int i = 0; i < diracVectorLength; ++i) {
			clover[0][i] = 0;
			clover[1][i] = 0;
			clover[2][i] = 0;
			clover[3][i] = 0;
			for (int j = 0; j < diracVectorLength; ++j) {
				clover[0][i] += I*FA[0][i][j]*input[site][0][j] + FB[0][i][j]*input[site][1][j];
				clover[1][i] += I*FA[1][i][j]*input[site][1][j] + FB[1][i][j]*input[site][0][j];
				clover[2][i] += I*FA[2][i][j]*input[site][2][j] + FB[2][i][j]*input[site][3][j];
				clover[3][i] += I*FA[3][i][j]*input[site][3][j] + FB[3][i][j]*input[site][2][j];
			}
		}

		if (gamma5) {
			//The final result is gamma5*input - kappa*hopping + kappa*csw*clover
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = + kappa*csw*clover[0][n] + std::complex<real_t>(real(input[site][0][n])-kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = + kappa*csw*clover[1][n] + std::complex<real_t>(real(input[site][1][n])-kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = + kappa*csw*clover[2][n] + std::complex<real_t>(-(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n]))),-(imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n]))));
				output[site][3][n] = + kappa*csw*clover[3][n] + std::complex<real_t>(-(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n]))),-(imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n]))));
			}
		}
		else {
			//The final result is input - kappa*gamma5*hopping + kappa*csw*gamma5*clover
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = + kappa*csw*clover[0][n] + std::complex<real_t>(real(input[site][0][n]) - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = + kappa*csw*clover[1][n] + std::complex<real_t>(real(input[site][1][n]) - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = - kappa*csw*clover[2][n] + std::complex<real_t>(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])),imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
				output[site][3][n] = - kappa*csw*clover[3][n] + std::complex<real_t>(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])),imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
			}
		}
*/
	}
	output.updateHalo();//TODO is needed?
}

void ImprovedDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & vector1, const reduced_dirac_vector_t & vector2, const complex& alpha) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;

	const reduced_fermion_lattice_t& linkconf = (lattice);

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		//the best
		std::complex<real_t> projection_spinor_minus[diracVectorLength], projection_spinor_plus[diracVectorLength], tmm, tmp;
		{
			{
				const size_t site_down = Vector::sdn(site,0);
				const size_t site_up = Vector::sup(site,0);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(vector1[site_down][0][n].real()+vector1[site_down][3][n].imag(),vector1[site_down][0][n].imag()-vector1[site_down][3][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(vector1[site_up][0][n].real()-vector1[site_up][3][n].imag(),vector1[site_up][0][n].imag()+vector1[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][0](i,n);
					}
					output[site][0][i] = alpha*vector2[site][0][i] + vector1[site][0][i] - kappa*(tmp+tmm);
					output[site][3][i] = alpha*vector2[site][3][i] - vector1[site][3][i] + kappa*std::complex<real_t>(tmm.imag() - tmp.imag(),tmp.real() - tmm.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(vector1[site_down][1][n].real()+vector1[site_down][2][n].imag(), vector1[site_down][1][n].imag()-vector1[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(vector1[site_up][1][n].real()-vector1[site_up][2][n].imag(), vector1[site_up][1][n].imag()+vector1[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][0](i,n);
					}

					output[site][1][i] = alpha*vector2[site][1][i] + vector1[site][1][i]- kappa*(tmp+tmm);
					output[site][2][i] = alpha*vector2[site][2][i] - vector1[site][2][i]+ kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real()-tmm.real());//(positI(tmp) +negatI(tmm));
				}
			}
			{
				const size_t site_down = Vector::sdn(site,1);
				const size_t site_up = Vector::sup(site,1);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = vector1[site_down][0][n] - (vector1[site_down][3][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = vector1[site_up][0][n] + (vector1[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][1](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*(tmm-tmp);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = vector1[site_down][1][n] + (vector1[site_down][2][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = vector1[site_up][1][n] - (vector1[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][1](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*(tmp-tmm);
				}
			}
			{
				const size_t site_down = Vector::sdn(site,2);
				const size_t site_up = Vector::sup(site,2);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(vector1[site_down][0][n].real() + vector1[site_down][2][n].imag(), vector1[site_down][0][n].imag() - vector1[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(vector1[site_up][0][n].real() - vector1[site_up][2][n].imag(), vector1[site_up][0][n].imag() + vector1[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][2](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real() - tmm.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = std::complex<real_t>(vector1[site_down][1][n].real() - vector1[site_down][3][n].imag(), vector1[site_down][1][n].imag() + vector1[site_down][3][n].real());
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = std::complex<real_t>(vector1[site_up][1][n].real() + vector1[site_up][3][n].imag(), vector1[site_up][1][n].imag() - vector1[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][2](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*std::complex<real_t>(tmp.imag() - tmm.imag(), tmm.real() - tmp.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,3);
				const size_t site_up = Vector::sup(site,3);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = vector1[site_down][0][n] + (vector1[site_down][2][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = vector1[site_up][0][n] - (vector1[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][3](i,n);
					}

					output[site][0][i] -= kappa*(tmp+tmm);
					output[site][2][i] += kappa*(tmp-tmm);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus[n] = vector1[site_down][1][n] + (vector1[site_down][3][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus[n] = vector1[site_up][1][n] - (vector1[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp = 0;
					tmm = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp += projection_spinor_minus[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm += projection_spinor_plus[n] * linkconf[site][3](i,n);
					}

					output[site][1][i] -= kappa*(tmp+tmm);
					output[site][3][i] += kappa*(tmp-tmm);
				}
			}

			//We store the result of the clover term in an intermediate vector
			GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				clover[0][i] = 0;
				clover[1][i] = 0;
				clover[2][i] = 0;
				clover[3][i] = 0;
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += multiply_by_I((-F[site][0].at(i,j)+F[site][5].at(i,j))*vector1[site][0][j]);
					clover[1][i] += multiply_by_I((+F[site][0].at(i,j)-F[site][5].at(i,j))*vector1[site][1][j]);
					clover[2][i] += multiply_by_I((+F[site][0].at(i,j)+F[site][5].at(i,j))*vector1[site][2][j]);
					clover[3][i] += multiply_by_I((-F[site][0].at(i,j)-F[site][5].at(i,j))*vector1[site][3][j]);
				}
			}

#ifdef ADJOINT
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += std::complex<real_t>(+F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*vector1[site][1][j];
					clover[1][i] += std::complex<real_t>(-F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*vector1[site][0][j];
					clover[2][i] += std::complex<real_t>(-F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*vector1[site][3][j];
					clover[3][i] += std::complex<real_t>(+F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*vector1[site][2][j];
				}
			}
#endif
#ifndef ADJOINT
//We store the result of the clover term in an intermediate vector
//GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += (+F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*vector1[site][1][j];
					clover[1][i] += (-F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*vector1[site][0][j];
					clover[2][i] += (-F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*vector1[site][3][j];
					clover[3][i] += (+F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*vector1[site][2][j];
				}
			}
#endif

			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] += (kappa*csw)*clover[0][n];
				output[site][1][n] += (kappa*csw)*clover[1][n];
				output[site][2][n] += (kappa*csw)*clover[2][n];
				output[site][3][n] += (kappa*csw)*clover[3][n];
			}
		}
		if (!gamma5) {
			for (int i = 0; i < diracVectorLength; ++i) {
				output[site][2][i] = -output[site][2][i]+static_cast<real_t>(2)*alpha*vector2[site][2][i];
				output[site][3][i] = -output[site][3][i]+static_cast<real_t>(2)*alpha*vector2[site][3][i];
			}
		}
/*

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
			projection_spinor[0][0][n] = std::complex<real_t>(real(vector1[site_sup_0][0][n])-imag(vector1[site_sup_0][3][n]),imag(vector1[site_sup_0][0][n])+real(vector1[site_sup_0][3][n]));
			projection_spinor[0][1][n] = std::complex<real_t>(real(vector1[site_sup_0][1][n])-imag(vector1[site_sup_0][2][n]),imag(vector1[site_sup_0][1][n])+real(vector1[site_sup_0][2][n]));
			projection_spinor[1][0][n] = std::complex<real_t>(real(vector1[site_sup_1][0][n])+real(vector1[site_sup_1][3][n]),imag(vector1[site_sup_1][0][n])+imag(vector1[site_sup_1][3][n]));
			projection_spinor[1][1][n] = std::complex<real_t>(real(vector1[site_sup_1][1][n])-real(vector1[site_sup_1][2][n]),imag(vector1[site_sup_1][1][n])-imag(vector1[site_sup_1][2][n]));
			projection_spinor[2][0][n] = std::complex<real_t>(real(vector1[site_sup_2][0][n])-imag(vector1[site_sup_2][2][n]),imag(vector1[site_sup_2][0][n])+real(vector1[site_sup_2][2][n]));
			projection_spinor[2][1][n] = std::complex<real_t>(real(vector1[site_sup_2][1][n])+imag(vector1[site_sup_2][3][n]),imag(vector1[site_sup_2][1][n])-real(vector1[site_sup_2][3][n]));
			projection_spinor[3][0][n] = std::complex<real_t>(real(vector1[site_sup_3][0][n])-real(vector1[site_sup_3][2][n]),imag(vector1[site_sup_3][0][n])-imag(vector1[site_sup_3][2][n]));
			projection_spinor[3][1][n] = std::complex<real_t>(real(vector1[site_sup_3][1][n])-real(vector1[site_sup_3][3][n]),imag(vector1[site_sup_3][1][n])-imag(vector1[site_sup_3][3][n]));
		}

		//Now we can put U(x,mu)*vector1(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 2; ++nu) {
				tmp_plus[mu][nu] = lattice[site][mu]*projection_spinor[mu][nu];
			}
		}

		int site_down_0 = Vector::sdn(site,0);
		int site_down_1 = Vector::sdn(site,1);
		int site_down_2 = Vector::sdn(site,2);
		int site_down_3 = Vector::sdn(site,3);

		for (int n = 0; n < diracVectorLength; ++n) {
			projection_spinor[0][0][n] = std::complex<real_t>(real(vector1[site_down_0][0][n])+imag(vector1[site_down_0][3][n]),imag(vector1[site_down_0][0][n])-real(vector1[site_down_0][3][n]));
			projection_spinor[0][1][n] = std::complex<real_t>(real(vector1[site_down_0][1][n])+imag(vector1[site_down_0][2][n]),imag(vector1[site_down_0][1][n])-real(vector1[site_down_0][2][n]));
			projection_spinor[1][0][n] = std::complex<real_t>(real(vector1[site_down_1][0][n])-real(vector1[site_down_1][3][n]),imag(vector1[site_down_1][0][n])-imag(vector1[site_down_1][3][n]));
			projection_spinor[1][1][n] = std::complex<real_t>(real(vector1[site_down_1][1][n])+real(vector1[site_down_1][2][n]),imag(vector1[site_down_1][1][n])+imag(vector1[site_down_1][2][n]));
			projection_spinor[2][0][n] = std::complex<real_t>(real(vector1[site_down_2][0][n])+imag(vector1[site_down_2][2][n]),imag(vector1[site_down_2][0][n])-real(vector1[site_down_2][2][n]));
			projection_spinor[2][1][n] = std::complex<real_t>(real(vector1[site_down_2][1][n])-imag(vector1[site_down_2][3][n]),imag(vector1[site_down_2][1][n])+real(vector1[site_down_2][3][n]));
			projection_spinor[3][0][n] = std::complex<real_t>(real(vector1[site_down_3][0][n])+real(vector1[site_down_3][2][n]),imag(vector1[site_down_3][0][n])+imag(vector1[site_down_3][2][n]));
			projection_spinor[3][1][n] = std::complex<real_t>(real(vector1[site_down_3][1][n])+real(vector1[site_down_3][3][n]),imag(vector1[site_down_3][1][n])+imag(vector1[site_down_3][3][n]));
		}

		//Then we put U(x-mu,mu)*vector1(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 2; ++nu) {
				GaugeVector tmp = htrans(lattice[Lattice::sdn(site,mu)][mu])*projection_spinor[mu][nu];
				tmp_minus[mu][nu] = tmp_plus[mu][nu] - tmp;
				tmp_plus[mu][nu] += tmp;
			}
		}

		//Clover matrix needed
		std::complex<real_t> FA[4][diracVectorLength][diracVectorLength];//TODO why real or complex?
		for (int i = 0; i < diracVectorLength; ++i) {
			for (int j = 0; j < diracVectorLength; ++j) {
				FA[0][i][j] = -F[site][0].at(i,j)+F[site][5].at(i,j);
				FA[1][i][j] = +F[site][0].at(i,j)-F[site][5].at(i,j);
				FA[2][i][j] = +F[site][0].at(i,j)+F[site][5].at(i,j);
				FA[3][i][j] = -F[site][0].at(i,j)-F[site][5].at(i,j);
			}
		}

		std::complex<real_t> FB[4][diracVectorLength][diracVectorLength];
		for (int i = 0; i < diracVectorLength; ++i) {
			for (int j = 0; j < diracVectorLength; ++j) {
				FB[0][i][j] = +F[site][1].at(i,j)+F[site][4].at(i,j)+I*(+F[site][2].at(i,j)-F[site][3].at(i,j));
				FB[1][i][j] = -F[site][1].at(i,j)-F[site][4].at(i,j)+I*(+F[site][2].at(i,j)-F[site][3].at(i,j));
				FB[2][i][j] = -F[site][1].at(i,j)+F[site][4].at(i,j)+I*(+F[site][2].at(i,j)+F[site][3].at(i,j));
				FB[3][i][j] = +F[site][1].at(i,j)-F[site][4].at(i,j)+I*(+F[site][2].at(i,j)+F[site][3].at(i,j));
			}
		}

		//We store the result of the clover term in an intermediate vector
		GaugeVector clover[4];
		for (int i = 0; i < diracVectorLength; ++i) {
			clover[0][i] = 0;
			clover[1][i] = 0;
			clover[2][i] = 0;
			clover[3][i] = 0;
			for (int j = 0; j < diracVectorLength; ++j) {
				clover[0][i] += I*FA[0][i][j]*vector1[site][0][j] + FB[0][i][j]*vector1[site][1][j];
				clover[1][i] += I*FA[1][i][j]*vector1[site][1][j] + FB[1][i][j]*vector1[site][0][j];
				clover[2][i] += I*FA[2][i][j]*vector1[site][2][j] + FB[2][i][j]*vector1[site][3][j];
				clover[3][i] += I*FA[3][i][j]*vector1[site][3][j] + FB[3][i][j]*vector1[site][2][j];
			}
		}

		if (gamma5) {
			//The final result is gamma5*input - kappa*hopping + kappa*csw*clover
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = alpha*vector2[site][0][n] + kappa*csw*clover[0][n] + std::complex<real_t>(real(vector1[site][0][n])-kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(vector1[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = alpha*vector2[site][1][n] + kappa*csw*clover[1][n] + std::complex<real_t>(real(vector1[site][1][n])-kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(vector1[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = alpha*vector2[site][2][n] + kappa*csw*clover[2][n] + std::complex<real_t>(-(real(vector1[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n]))),-(imag(vector1[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n]))));
				output[site][3][n] = alpha*vector2[site][3][n] + kappa*csw*clover[3][n] + std::complex<real_t>(-(real(vector1[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n]))),-(imag(vector1[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n]))));
			}
		}
		else {
			//The final result is input - kappa*gamma5*hopping + kappa*csw*gamma5*clover
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = alpha*vector2[site][0][n] + kappa*csw*clover[0][n] + std::complex<real_t>(real(vector1[site][0][n]) - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(vector1[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = alpha*vector2[site][1][n] + kappa*csw*clover[1][n] + std::complex<real_t>(real(vector1[site][1][n]) - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(vector1[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = alpha*vector2[site][2][n] - kappa*csw*clover[2][n] + std::complex<real_t>(real(vector1[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])),imag(vector1[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
				output[site][3][n] = alpha*vector2[site][3][n] - kappa*csw*clover[3][n] + std::complex<real_t>(real(vector1[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])),imag(vector1[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
			}
		}
		*/
	}
	output.updateHalo();//TODO is needed?
}

void ImprovedDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	this->updateFieldStrength(_lattice);
}

real_t ImprovedDiracWilsonOperator::getCSW() const {
	return csw;
}

void ImprovedDiracWilsonOperator::setCSW(real_t _csw) {
	csw = _csw;
}

void ImprovedDiracWilsonOperator::updateFieldStrength(const extended_fermion_lattice_t& _lattice) {
	typedef extended_fermion_lattice_t LT;
	extended_field_strength_lattice_t tmpF;
	//We must calculate everything in extended layout!!!
#pragma omp parallel for
	for (int site = 0; site < _lattice.localsize; ++site) {
		tmpF[site][0] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 1)][1])*(_lattice[LT::sdn(LT::sdn(site, 0), 1)][0])*(_lattice[LT::sdn(site, 1)][1]) + htrans(_lattice[LT::sdn(site, 1)][1])*(_lattice[LT::sdn(site, 1)][0])*(_lattice[LT::sup(LT::sdn(site, 1), 0)][1])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][1])*htrans(_lattice[LT::sup(site, 1)][0])*htrans(_lattice[site][1]) + (_lattice[site][1])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 1)][0])*htrans(_lattice[LT::sdn(site, 0)][1])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[site][1] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 2)][2])*(_lattice[LT::sdn(LT::sdn(site, 0), 2)][0])*(_lattice[LT::sdn(site, 2)][2]) + htrans(_lattice[LT::sdn(site, 2)][2])*(_lattice[LT::sdn(site, 2)][0])*(_lattice[LT::sup(LT::sdn(site, 2), 0)][2])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][2])*htrans(_lattice[LT::sup(site, 2)][0])*htrans(_lattice[site][2]) + (_lattice[site][2])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 2)][0])*htrans(_lattice[LT::sdn(site, 0)][2])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[site][2] = htrans(_lattice[LT::sdn(site, 0)][0])*htrans(_lattice[LT::sdn(LT::sdn(site, 0), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 0), 3)][0])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][0])*(_lattice[LT::sup(LT::sdn(site, 3), 0)][3])*htrans(_lattice[site][0]) + (_lattice[site][0])*(_lattice[LT::sup(site, 0)][3])*htrans(_lattice[LT::sup(site, 3)][0])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 0), 3)][0])*htrans(_lattice[LT::sdn(site, 0)][3])*(_lattice[LT::sdn(site, 0)][0]);
		tmpF[site][3] = htrans(_lattice[LT::sdn(site, 1)][1])*htrans(_lattice[LT::sdn(LT::sdn(site, 1), 2)][2])*(_lattice[LT::sdn(LT::sdn(site, 1), 2)][1])*(_lattice[LT::sdn(site, 2)][2]) + htrans(_lattice[LT::sdn(site, 2)][2])*(_lattice[LT::sdn(site, 2)][1])*(_lattice[LT::sup(LT::sdn(site, 2), 1)][2])*htrans(_lattice[site][1]) + (_lattice[site][1])*(_lattice[LT::sup(site, 1)][2])*htrans(_lattice[LT::sup(site, 2)][1])*htrans(_lattice[site][2]) + (_lattice[site][2])*htrans(_lattice[LT::sup(LT::sdn(site, 1), 2)][1])*htrans(_lattice[LT::sdn(site, 1)][2])*(_lattice[LT::sdn(site, 1)][1]);
		tmpF[site][4] = htrans(_lattice[LT::sdn(site, 1)][1])*htrans(_lattice[LT::sdn(LT::sdn(site, 1), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 1), 3)][1])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][1])*(_lattice[LT::sup(LT::sdn(site, 3), 1)][3])*htrans(_lattice[site][1]) + (_lattice[site][1])*(_lattice[LT::sup(site, 1)][3])*htrans(_lattice[LT::sup(site, 3)][1])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 1), 3)][1])*htrans(_lattice[LT::sdn(site, 1)][3])*(_lattice[LT::sdn(site, 1)][1]);
		tmpF[site][5] = htrans(_lattice[LT::sdn(site, 2)][2])*htrans(_lattice[LT::sdn(LT::sdn(site, 2), 3)][3])*(_lattice[LT::sdn(LT::sdn(site, 2), 3)][2])*(_lattice[LT::sdn(site, 3)][3]) + htrans(_lattice[LT::sdn(site, 3)][3])*(_lattice[LT::sdn(site, 3)][2])*(_lattice[LT::sup(LT::sdn(site, 3), 2)][3])*htrans(_lattice[site][2]) + (_lattice[site][2])*(_lattice[LT::sup(site, 2)][3])*htrans(_lattice[LT::sup(site, 3)][2])*htrans(_lattice[site][3]) + (_lattice[site][3])*htrans(_lattice[LT::sup(LT::sdn(site, 2), 3)][2])*htrans(_lattice[LT::sdn(site, 2)][3])*(_lattice[LT::sdn(site, 2)][2]);
		for (unsigned int i = 0; i < 6; ++i) {
			//Manual antialiasing, error of eigen!
			FermionicGroup antialias = tmpF[site][i];
			tmpF[site][i] = (1./8.)*(antialias - htrans(antialias));
		}
	}
	tmpF.updateHalo();
	//Now we get the reduced lattice
	F = tmpF;
}

FermionForce* ImprovedDiracWilsonOperator::getForce() const {
	return new ImprovedFermionForce(kappa, csw);
}

} /* namespace Update */
