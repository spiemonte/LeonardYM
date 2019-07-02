/*
 * ImprovedDiracWilsonOperator.cpp
 *
 *  Created on: May 4, 2012
 *      Author: spiem_01
 */

#include "ImprovedDiracWilsonOperator.h"
#include "hmc_forces/ImprovedFermionForce.h"

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
					output[site][2][i] = -input[site][2][i]+ kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real()-tmm.real());
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
	}
	output.updateHalo();
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
					output[site][2][i] = alpha*vector2[site][2][i] - vector1[site][2][i]+ kappa*std::complex<real_t>(tmm.imag()-tmp.imag(),tmp.real()-tmm.real());
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
	}
	output.updateHalo();
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
