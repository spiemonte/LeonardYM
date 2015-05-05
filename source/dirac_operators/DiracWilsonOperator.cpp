/*
 * DiracWilsonOperator.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: spiem_01
 */

#include "DiracWilsonOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"
/*#include <x86intrin.h>
#include <immintrin.h>
#include <malloc.h>*/

namespace Update {

DiracWilsonOperator::DiracWilsonOperator() : DiracOperator() { }

DiracWilsonOperator::DiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

DiracWilsonOperator::~DiracWilsonOperator() { }

/*void multiply(double* output, std::complex<double>* input, double* lattice, double& kappa) {
	typedef adjoint_lattice_t LT;
	double* llRe = (double*) memalign(64,sizeof(double)*4);
	double* llIm = (double*) memalign(64,sizeof(double)*4);
	for (unsigned int site = 0; site < LT::localsize; ++site) {
		std::complex<double> tmp[4][4][3];

		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_up = LT::sup(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				__m256d row1 = _mm256_load_pd(lattice[9*(4*site+mu)]);
				__m256d row2 = _mm256_load_pd(lattice[9*(4*site+mu)+3]);
				__m256d row3 = _mm256_load_pd(lattice[9*(4*site+mu)+6]);

				//Then we perform the multiplication
				__m256d res1Re = _mm256_mul_pd(row1,_mm256_set1_pd(real(input[site_up][nu][0])));
				__m256d res2Re = _mm256_mul_pd(row2,_mm256_set1_pd(real(input[site_up][nu][1])));
				__m256d res3Re = _mm256_mul_pd(row3,_mm256_set1_pd(real(input[site_up][nu][2])));
				__m256d res1Im = _mm256_mul_pd(row1,_mm256_set1_pd(imag(input[site_up][nu][0])));
				__m256d res2Im = _mm256_mul_pd(row2,_mm256_set1_pd(imag(input[site_up][nu][1])));
				__m256d res3Im = _mm256_mul_pd(row3,_mm256_set1_pd(imag(input[site_up][nu][2])));

				//Now we get the final result
				__m256d resfRe = _mm256_add_pd(_mm256_add_pd(res1Re,res2Re),res3Re);
				__m256d resfIm = _mm256_add_pd(_mm256_add_pd(res1Im,res2Im),res3Im);

				//We get the result
				_mm256_store_pd(llRe,resfRe);
				_mm256_store_pd(llIm,resfIm);
				//Everything is stored in the reverse order!
				tmp[mu][nu][0] = std::complex<real_t>(llRe[3],llIm[3]);
				tmp[mu][nu][1] = std::complex<real_t>(llRe[2],llIm[2]);
				tmp[mu][nu][3] = std::complex<real_t>(llRe[1],llIm[1]);
			}
		}
	}
}*/

/*void DiracWilsonOperator::multiply(dirac_vector_t& output, const dirac_vector_t& input) {
	typedef adjoint_lattice_t LT;
	typedef dirac_vector_t DV;
#ifdef MULTITHREADING
	#pragma omp parallel for
#endif
	for (unsigned int site = 0; site < lattice->localsize; ++site) {
		//First we start the hopping parameter terms
		GaugeVector tmp[4][4];

		//First we put U(x,mu)*input(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_up_vector = DV::sup(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = (*lattice)[site][mu]*input[site_up_vector][nu];
			}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4];
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] = (tmp[0][0][n]+I*tmp[0][3][n]+tmp[1][0][n]+tmp[1][3][n]+tmp[2][0][n]+I*tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[1][n] = (tmp[0][1][n]+I*tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+tmp[2][1][n]-I*tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
			hopping[2][n] = (I*tmp[0][1][n]-tmp[0][2][n]+tmp[1][1][n]-tmp[1][2][n]+I*tmp[2][0][n]-tmp[2][2][n]+tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] = (I*tmp[0][0][n]-tmp[0][3][n]-tmp[1][0][n]-tmp[1][3][n]-I*tmp[2][1][n]-tmp[2][3][n]+tmp[3][1][n]-tmp[3][3][n]);
		}

		//Then we put U(x-mu,mu)*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_down_lattice = LT::sdn(site,mu);
			int site_down_vector = DV::sdn(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = htrans((*lattice)[site_down_lattice][mu])*input[site_down_vector][nu];
			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			hopping[0][n] += (tmp[0][0][n]-I*tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+tmp[2][0][n]-I*tmp[2][2][n]+tmp[3][0][n]+tmp[3][2][n]);
			hopping[1][n] += (tmp[0][1][n]-I*tmp[0][2][n]+tmp[1][1][n]+tmp[1][2][n]+tmp[2][1][n]+I*tmp[2][3][n]+tmp[3][1][n]+tmp[3][3][n]);
			hopping[2][n] += (-I*tmp[0][1][n]-tmp[0][2][n]-tmp[1][1][n]-tmp[1][2][n]-I*tmp[2][0][n]-tmp[2][2][n]-tmp[3][0][n]-tmp[3][2][n]);
			hopping[3][n] += (-I*tmp[0][0][n]-tmp[0][3][n]+tmp[1][0][n]-tmp[1][3][n]+I*tmp[2][1][n]-tmp[2][3][n]-tmp[3][1][n]-tmp[3][3][n]);
		}

		//The final result is gamma5*input - kappa*hopping
		for (unsigned int n = 0; n < diracVectorLength; ++n) {
			output[site][0][n] = input[site][0][n] - kappa*hopping[0][n];
			output[site][1][n] = input[site][1][n] - kappa*hopping[1][n];
			output[site][2][n] = -(input[site][2][n] + kappa*hopping[2][n]);
			output[site][3][n] = -(input[site][3][n] + kappa*hopping[3][n]);
		}
	}
	output.updateHalo();
}*/

inline GaugeVector multiply_by_I(const GaugeVector& a) {
	GaugeVector result;
	for (unsigned int i = 0; i < diracVectorLength; ++i) result[i] = std::complex<real_t>(-imag(a[i]),real(a[i]));
	return result;
}

inline real_t conj(const real_t& t) {
	return t;
}

void DiracWilsonOperator::multiply(reduced_dirac_vector_t& output1, reduced_dirac_vector_t& output2, reduced_dirac_vector_t& output3, reduced_dirac_vector_t& output4, const reduced_dirac_vector_t& input1, const reduced_dirac_vector_t& input2, const reduced_dirac_vector_t& input3, const reduced_dirac_vector_t& input4) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;
	const reduced_fermion_lattice_t& linkconf = (lattice);

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		/*std::complex<real_t> projection_spinor_minus1[diracVectorLength], projection_spinor_plus1[diracVectorLength], tmm1, tmp1;
		std::complex<real_t> projection_spinor_minus2[diracVectorLength], projection_spinor_plus2[diracVectorLength], tmm2, tmp2;
		std::complex<real_t> projection_spinor_minus3[diracVectorLength], projection_spinor_plus3[diracVectorLength], tmm3, tmp3;
		std::complex<real_t> projection_spinor_minus4[diracVectorLength], projection_spinor_plus4[diracVectorLength], tmm4, tmp4;
		{
			{
				const size_t site_down = Vector::sdn(site,0);
				const size_t site_up = Vector::sup(site,0);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real()+input1[site_down][3][n].imag(),input1[site_down][0][n].imag()-input1[site_down][3][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real()+input2[site_down][3][n].imag(),input2[site_down][0][n].imag()-input2[site_down][3][n].real());
					projection_spinor_minus3[n] = std::complex<real_t>(input3[site_down][0][n].real()+input3[site_down][3][n].imag(),input3[site_down][0][n].imag()-input3[site_down][3][n].real());
					projection_spinor_minus4[n] = std::complex<real_t>(input4[site_down][0][n].real()+input4[site_down][3][n].imag(),input4[site_down][0][n].imag()-input4[site_down][3][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real()-input1[site_up][3][n].imag(),input1[site_up][0][n].imag()+input1[site_up][3][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real()-input2[site_up][3][n].imag(),input2[site_up][0][n].imag()+input2[site_up][3][n].real());
					projection_spinor_plus3[n] = std::complex<real_t>(input3[site_up][0][n].real()-input3[site_up][3][n].imag(),input3[site_up][0][n].imag()+input3[site_up][3][n].real());
					projection_spinor_plus4[n] = std::complex<real_t>(input4[site_up][0][n].real()-input4[site_up][3][n].imag(),input4[site_up][0][n].imag()+input4[site_up][3][n].real());

				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));

					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][0](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][0](i,n);
					}
					output1[site][0][i] = input1[site][0][i] - kappa*(tmp1+tmm1);
					output2[site][0][i] = input2[site][0][i] - kappa*(tmp2+tmm2);
					output3[site][0][i] = input3[site][0][i] - kappa*(tmp3+tmm3);
					output4[site][0][i] = input4[site][0][i] - kappa*(tmp4+tmm4);

					output1[site][3][i] = -input1[site][3][i] + kappa*std::complex<real_t>(tmm1.imag() - tmp1.imag(),tmp1.real() - tmm1.real());
					output2[site][3][i] = -input2[site][3][i] + kappa*std::complex<real_t>(tmm2.imag() - tmp2.imag(),tmp2.real() - tmm2.real());
					output3[site][3][i] = -input3[site][3][i] + kappa*std::complex<real_t>(tmm3.imag() - tmp3.imag(),tmp3.real() - tmm3.real());
					output4[site][3][i] = -input4[site][3][i] + kappa*std::complex<real_t>(tmm4.imag() - tmp4.imag(),tmp4.real() - tmm4.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real()+input1[site_down][2][n].imag(), input1[site_down][1][n].imag()-input1[site_down][2][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real()+input2[site_down][2][n].imag(), input2[site_down][1][n].imag()-input2[site_down][2][n].real());
					projection_spinor_minus3[n] = std::complex<real_t>(input3[site_down][1][n].real()+input3[site_down][2][n].imag(), input3[site_down][1][n].imag()-input3[site_down][2][n].real());
					projection_spinor_minus4[n] = std::complex<real_t>(input4[site_down][1][n].real()+input4[site_down][2][n].imag(), input4[site_down][1][n].imag()-input4[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real()-input1[site_up][2][n].imag(),input1[site_up][1][n].imag()+input1[site_up][2][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real()-input2[site_up][2][n].imag(),input2[site_up][1][n].imag()+input2[site_up][2][n].real());
					projection_spinor_plus3[n] = std::complex<real_t>(input3[site_up][1][n].real()-input3[site_up][2][n].imag(),input3[site_up][1][n].imag()+input3[site_up][2][n].real());
					projection_spinor_plus4[n] = std::complex<real_t>(input4[site_up][1][n].real()-input4[site_up][2][n].imag(),input4[site_up][1][n].imag()+input4[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));

					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][0](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][0](i,n);
					}

					output1[site][1][i] = input1[site][1][i]- kappa*(tmp1+tmm1);
					output2[site][1][i] = input2[site][1][i]- kappa*(tmp2+tmm2);
					output3[site][1][i] = input3[site][1][i]- kappa*(tmp3+tmm3);
					output4[site][1][i] = input4[site][1][i]- kappa*(tmp4+tmm4);
					output1[site][2][i] = -input1[site][2][i]+ kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real()-tmm1.real());
					output2[site][2][i] = -input2[site][2][i]+ kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real()-tmm2.real());
					output3[site][2][i] = -input3[site][2][i]+ kappa*std::complex<real_t>(tmm3.imag()-tmp3.imag(),tmp3.real()-tmm3.real());
					output4[site][2][i] = -input4[site][2][i]+ kappa*std::complex<real_t>(tmm4.imag()-tmp4.imag(),tmp4.real()-tmm4.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,1);
				const size_t site_up = Vector::sup(site,1);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][0][n] - (input1[site_down][3][n]);
					projection_spinor_minus2[n] = input2[site_down][0][n] - (input2[site_down][3][n]);
					projection_spinor_minus3[n] = input3[site_down][0][n] - (input3[site_down][3][n]);
					projection_spinor_minus4[n] = input4[site_down][0][n] - (input4[site_down][3][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][0][n] + (input1[site_up][3][n]);
					projection_spinor_plus2[n] = input2[site_up][0][n] + (input2[site_up][3][n]);
					projection_spinor_plus3[n] = input3[site_up][0][n] + (input3[site_up][3][n]);
					projection_spinor_plus4[n] = input4[site_up][0][n] + (input4[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][1](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][1](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output3[site][0][i] -= kappa*(tmp3+tmm3);
					output4[site][0][i] -= kappa*(tmp4+tmm4);

					output1[site][3][i] += kappa*(tmm1-tmp1);
					output2[site][3][i] += kappa*(tmm2-tmp2);
					output3[site][3][i] += kappa*(tmm3-tmp3);
					output4[site][3][i] += kappa*(tmm4-tmp4);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][2][n]);
					projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][2][n]);
					projection_spinor_minus3[n] = input3[site_down][1][n] + (input3[site_down][2][n]);
					projection_spinor_minus4[n] = input4[site_down][1][n] + (input4[site_down][2][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input2[site_up][1][n] - (input1[site_up][2][n]);
					projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][2][n]);
					projection_spinor_plus3[n] = input3[site_up][1][n] - (input3[site_up][2][n]);
					projection_spinor_plus4[n] = input4[site_up][1][n] - (input4[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][1](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][1](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output3[site][1][i] -= kappa*(tmp3+tmm3);
					output4[site][1][i] -= kappa*(tmp4+tmm4);
					output1[site][2][i] += kappa*(tmp1-tmm1);
					output2[site][2][i] += kappa*(tmp2-tmm2);
					output3[site][2][i] += kappa*(tmp3-tmm3);
					output4[site][2][i] += kappa*(tmp4-tmm4);
				}
			}
			{
				const size_t site_down = Vector::sdn(site,2);
				const size_t site_up = Vector::sup(site,2);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real() + input1[site_down][2][n].imag(), input1[site_down][0][n].imag() - input1[site_down][2][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real() + input2[site_down][2][n].imag(), input2[site_down][0][n].imag() - input2[site_down][2][n].real());
					projection_spinor_minus3[n] = std::complex<real_t>(input3[site_down][0][n].real() + input3[site_down][2][n].imag(), input3[site_down][0][n].imag() - input3[site_down][2][n].real());
					projection_spinor_minus4[n] = std::complex<real_t>(input4[site_down][0][n].real() + input4[site_down][2][n].imag(), input4[site_down][0][n].imag() - input4[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real() - input1[site_up][2][n].imag(), input1[site_up][0][n].imag() + input1[site_up][2][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real() - input2[site_up][2][n].imag(), input2[site_up][0][n].imag() + input2[site_up][2][n].real());
					projection_spinor_plus3[n] = std::complex<real_t>(input3[site_up][0][n].real() - input3[site_up][2][n].imag(), input3[site_up][0][n].imag() + input3[site_up][2][n].real());
					projection_spinor_plus4[n] = std::complex<real_t>(input4[site_up][0][n].real() - input4[site_up][2][n].imag(), input4[site_up][0][n].imag() + input4[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][2](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][2](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output3[site][0][i] -= kappa*(tmp3+tmm3);
					output4[site][0][i] -= kappa*(tmp4+tmm4);
					output1[site][2][i] += kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real() - tmm1.real());
					output2[site][2][i] += kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real() - tmm2.real());
					output3[site][2][i] += kappa*std::complex<real_t>(tmm3.imag()-tmp3.imag(),tmp3.real() - tmm3.real());
					output4[site][2][i] += kappa*std::complex<real_t>(tmm4.imag()-tmp4.imag(),tmp4.real() - tmm4.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real() - input1[site_down][3][n].imag(), input1[site_down][1][n].imag() + input1[site_down][3][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real() - input2[site_down][3][n].imag(), input2[site_down][1][n].imag() + input2[site_down][3][n].real());
					projection_spinor_minus3[n] = std::complex<real_t>(input3[site_down][1][n].real() - input3[site_down][3][n].imag(), input3[site_down][1][n].imag() + input3[site_down][3][n].real());
					projection_spinor_minus4[n] = std::complex<real_t>(input4[site_down][1][n].real() - input4[site_down][3][n].imag(), input4[site_down][1][n].imag() + input4[site_down][3][n].real());
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real() + input1[site_up][3][n].imag(), input1[site_up][1][n].imag() - input1[site_up][3][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real() + input2[site_up][3][n].imag(), input2[site_up][1][n].imag() - input2[site_up][3][n].real());
					projection_spinor_plus3[n] = std::complex<real_t>(input3[site_up][1][n].real() + input3[site_up][3][n].imag(), input3[site_up][1][n].imag() - input3[site_up][3][n].real());
					projection_spinor_plus4[n] = std::complex<real_t>(input4[site_up][1][n].real() + input4[site_up][3][n].imag(), input4[site_up][1][n].imag() - input4[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][2](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][2](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output3[site][1][i] -= kappa*(tmp3+tmm3);
					output4[site][1][i] -= kappa*(tmp4+tmm4);
					output1[site][3][i] += kappa*std::complex<real_t>(tmp1.imag() - tmm1.imag(), tmm1.real() - tmp1.real());
					output2[site][3][i] += kappa*std::complex<real_t>(tmp2.imag() - tmm2.imag(), tmm2.real() - tmp2.real());
					output3[site][3][i] += kappa*std::complex<real_t>(tmp3.imag() - tmm3.imag(), tmm3.real() - tmp3.real());
					output4[site][3][i] += kappa*std::complex<real_t>(tmp4.imag() - tmm4.imag(), tmm4.real() - tmp4.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,3);
				const size_t site_up = Vector::sup(site,3);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][0][n] + (input1[site_down][2][n]);
					projection_spinor_minus2[n] = input2[site_down][0][n] + (input2[site_down][2][n]);
					projection_spinor_minus3[n] = input3[site_down][0][n] + (input3[site_down][2][n]);
					projection_spinor_minus4[n] = input4[site_down][0][n] + (input4[site_down][2][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][0][n] - (input1[site_up][2][n]);
					projection_spinor_plus2[n] = input2[site_up][0][n] - (input2[site_up][2][n]);
					projection_spinor_plus3[n] = input3[site_up][0][n] - (input3[site_up][2][n]);
					projection_spinor_plus4[n] = input4[site_up][0][n] - (input4[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][3](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][3](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output3[site][0][i] -= kappa*(tmp3+tmm3);
					output4[site][0][i] -= kappa*(tmp4+tmm4);
					output1[site][2][i] += kappa*(tmp1-tmm1);
					output2[site][2][i] += kappa*(tmp2-tmm2);
					output3[site][2][i] += kappa*(tmp3-tmm3);
					output4[site][2][i] += kappa*(tmp4-tmm4);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][3][n]);
					projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][3][n]);
					projection_spinor_minus3[n] = input3[site_down][1][n] + (input3[site_down][3][n]);
					projection_spinor_minus4[n] = input4[site_down][1][n] + (input4[site_down][3][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][1][n] - (input1[site_up][3][n]);
					projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][3][n]);
					projection_spinor_plus3[n] = input3[site_up][1][n] - (input3[site_up][3][n]);
					projection_spinor_plus4[n] = input4[site_up][1][n] - (input4[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmp3 = 0;
					tmp4 = 0;
					tmm1 = 0;
					tmm2 = 0;
					tmm3 = 0;
					tmm4 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp3 += projection_spinor_minus3[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp4 += projection_spinor_minus4[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						tmm3 += projection_spinor_plus3[n] * linkconf[site][3](i,n);
						tmm4 += projection_spinor_plus4[n] * linkconf[site][3](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output3[site][1][i] -= kappa*(tmp3+tmm3);
					output4[site][1][i] -= kappa*(tmp4+tmm4);
					output1[site][3][i] += kappa*(tmp1-tmm1);
					output2[site][3][i] += kappa*(tmp2-tmm2);
					output3[site][3][i] += kappa*(tmp3-tmm3);
					output4[site][3][i] += kappa*(tmp4-tmm4);
				}
			}
		}*/

		std::complex<real_t> projection_spinor_minus1[diracVectorLength], projection_spinor_plus1[diracVectorLength], tmm1, tmp1;
		std::complex<real_t> projection_spinor_minus2[diracVectorLength], projection_spinor_plus2[diracVectorLength], tmm2, tmp2;
		{
			{
				const size_t site_down = Vector::sdn(site,0);
				const size_t site_up = Vector::sup(site,0);
				{
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real()+input1[site_down][3][n].imag(),input1[site_down][0][n].imag()-input1[site_down][3][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real()+input2[site_down][3][n].imag(),input2[site_down][0][n].imag()-input2[site_down][3][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real()-input1[site_up][3][n].imag(),input1[site_up][0][n].imag()+input1[site_up][3][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real()-input2[site_up][3][n].imag(),input2[site_up][0][n].imag()+input2[site_up][3][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp1 = 0;
						tmm2 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						}
						output1[site][0][i] = input1[site][0][i] - kappa*(tmp1+tmm2);
						output2[site][0][i] = input2[site][0][i] - kappa*(tmp2+tmm2);
						output1[site][3][i] = -input1[site][3][i] + kappa*std::complex<real_t>(tmm1.imag() - tmp1.imag(),tmp1.real() - tmm1.real());
						output2[site][3][i] = -input2[site][3][i] + kappa*std::complex<real_t>(tmm2.imag() - tmp2.imag(),tmp2.real() - tmm2.real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real()+input1[site_down][2][n].imag(), input1[site_down][1][n].imag()-input1[site_down][2][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real()+input2[site_down][2][n].imag(), input2[site_down][1][n].imag()-input2[site_down][2][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real()-input1[site_up][2][n].imag(),input1[site_up][1][n].imag()+input1[site_up][2][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real()-input2[site_up][2][n].imag(),input2[site_up][1][n].imag()+input2[site_up][2][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						}

						output1[site][1][i] = input1[site][1][i]- kappa*(tmp1+tmm1);
						output2[site][1][i] = input2[site][1][i]- kappa*(tmp2+tmm2);
						output1[site][2][i] = -input1[site][2][i]+ kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real()-tmm1.real());
						output2[site][2][i] = -input2[site][2][i]+ kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real()-tmm2.real());
					}
				}
				{
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input3[site_down][0][n].real()+input3[site_down][3][n].imag(),input3[site_down][0][n].imag()-input3[site_down][3][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input4[site_down][0][n].real()+input4[site_down][3][n].imag(),input4[site_down][0][n].imag()-input4[site_down][3][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input3[site_up][0][n].real()-input3[site_up][3][n].imag(),input3[site_up][0][n].imag()+input3[site_up][3][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input4[site_up][0][n].real()-input4[site_up][3][n].imag(),input4[site_up][0][n].imag()+input4[site_up][3][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp1 = 0;
						tmm2 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						}
						output3[site][0][i] = input3[site][0][i] - kappa*(tmp1+tmm2);
						output4[site][0][i] = input4[site][0][i] - kappa*(tmp2+tmm2);
						output3[site][3][i] = -input3[site][3][i] + kappa*std::complex<real_t>(tmm1.imag() - tmp1.imag(),tmp1.real() - tmm1.real());
						output4[site][3][i] = -input4[site][3][i] + kappa*std::complex<real_t>(tmm2.imag() - tmp2.imag(),tmp2.real() - tmm2.real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input3[site_down][1][n].real()+input3[site_down][2][n].imag(), input3[site_down][1][n].imag()-input3[site_down][2][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input4[site_down][1][n].real()+input4[site_down][2][n].imag(), input4[site_down][1][n].imag()-input4[site_down][2][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input3[site_up][1][n].real()-input3[site_up][2][n].imag(),input3[site_up][1][n].imag()+input3[site_up][2][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input4[site_up][1][n].real()-input4[site_up][2][n].imag(),input4[site_up][1][n].imag()+input4[site_up][2][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
						}

						output3[site][1][i] = input3[site][1][i]- kappa*(tmp1+tmm1);
						output4[site][1][i] = input4[site][1][i]- kappa*(tmp2+tmm2);
						output3[site][2][i] = -input3[site][2][i]+ kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real()-tmm1.real());
						output4[site][2][i] = -input4[site][2][i]+ kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real()-tmm2.real());
					}
				}
			}




			{
				const size_t site_down = Vector::sdn(site,1);
				const size_t site_up = Vector::sup(site,1);
				{

					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input1[site_down][0][n] - (input1[site_down][3][n]);
						projection_spinor_minus2[n] = input2[site_down][0][n] - (input2[site_down][3][n]);
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input1[site_up][0][n] + (input1[site_up][3][n]);
						projection_spinor_plus2[n] = input2[site_up][0][n] + (input2[site_up][3][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for(int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						}
						for(int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						}

						output1[site][0][i] -= kappa*(tmp1+tmm1);
						output2[site][0][i] -= kappa*(tmp2+tmm2);
						output1[site][3][i] += kappa*(tmm1-tmp1);
						output2[site][3][i] += kappa*(tmm2-tmp2);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][2][n]);
						projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][2][n]);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input2[site_up][1][n] - (input1[site_up][2][n]);
						projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][2][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm2 = 0;
						tmm1 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						}

						output1[site][1][i] -= kappa*(tmp1+tmm1);
						output2[site][1][i] -= kappa*(tmp2+tmm2);
						output1[site][2][i] += kappa*(tmp1-tmm1);
						output2[site][2][i] += kappa*(tmp2-tmm2);
					}
				}
				{
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input3[site_down][0][n] - (input3[site_down][3][n]);
						projection_spinor_minus2[n] = input4[site_down][0][n] - (input4[site_down][3][n]);
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input3[site_up][0][n] + (input3[site_up][3][n]);
						projection_spinor_plus2[n] = input4[site_up][0][n] + (input4[site_up][3][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for(int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						}
						for(int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						}

						output3[site][0][i] -= kappa*(tmp1+tmm1);
						output4[site][0][i] -= kappa*(tmp2+tmm2);
						output3[site][3][i] += kappa*(tmm1-tmp1);
						output4[site][3][i] += kappa*(tmm2-tmp2);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input3[site_down][1][n] + (input3[site_down][2][n]);
						projection_spinor_minus2[n] = input4[site_down][1][n] + (input4[site_down][2][n]);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input4[site_up][1][n] - (input3[site_up][2][n]);
						projection_spinor_plus2[n] = input4[site_up][1][n] - (input4[site_up][2][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm2 = 0;
						tmm1 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
						}

						output3[site][1][i] -= kappa*(tmp1+tmm1);
						output4[site][1][i] -= kappa*(tmp2+tmm2);
						output3[site][2][i] += kappa*(tmp1-tmm1);
						output4[site][2][i] += kappa*(tmp2-tmm2);
					}
				}
			}




			{
				const size_t site_down = Vector::sdn(site,2);
				const size_t site_up = Vector::sup(site,2);
				{

					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real() + input1[site_down][2][n].imag(), input1[site_down][0][n].imag() - input1[site_down][2][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real() + input2[site_down][2][n].imag(), input2[site_down][0][n].imag() - input2[site_down][2][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real() - input1[site_up][2][n].imag(), input1[site_up][0][n].imag() + input1[site_up][2][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real() - input2[site_up][2][n].imag(), input2[site_up][0][n].imag() + input2[site_up][2][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						}

						output1[site][0][i] -= kappa*(tmp1+tmm1);
						output2[site][0][i] -= kappa*(tmp2+tmm2);
						output1[site][2][i] += kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real() - tmm1.real());
						output2[site][2][i] += kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real() - tmm2.real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real() - input1[site_down][3][n].imag(), input1[site_down][1][n].imag() + input1[site_down][3][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real() - input2[site_down][3][n].imag(), input2[site_down][1][n].imag() + input2[site_down][3][n].real());
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real() + input1[site_up][3][n].imag(), input1[site_up][1][n].imag() - input1[site_up][3][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real() + input2[site_up][3][n].imag(), input2[site_up][1][n].imag() - input2[site_up][3][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for(int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						}
						for(int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						}

						output1[site][1][i] -= kappa*(tmp1+tmm1);
						output2[site][1][i] -= kappa*(tmp2+tmm2);
						output1[site][3][i] += kappa*std::complex<real_t>(tmp1.imag() - tmm1.imag(), tmm1.real() - tmp1.real());
						output2[site][3][i] += kappa*std::complex<real_t>(tmp2.imag() - tmm2.imag(), tmm2.real() - tmp2.real());
					}
				}
				{
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input3[site_down][0][n].real() + input3[site_down][2][n].imag(), input3[site_down][0][n].imag() - input3[site_down][2][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input4[site_down][0][n].real() + input4[site_down][2][n].imag(), input4[site_down][0][n].imag() - input4[site_down][2][n].real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input3[site_up][0][n].real() - input3[site_up][2][n].imag(), input3[site_up][0][n].imag() + input3[site_up][2][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input4[site_up][0][n].real() - input4[site_up][2][n].imag(), input4[site_up][0][n].imag() + input4[site_up][2][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						}

						output3[site][0][i] -= kappa*(tmp1+tmm1);
						output4[site][0][i] -= kappa*(tmp2+tmm2);
						output3[site][2][i] += kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real() - tmm1.real());
						output4[site][2][i] += kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real() - tmm2.real());
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = std::complex<real_t>(input3[site_down][1][n].real() - input3[site_down][3][n].imag(), input3[site_down][1][n].imag() + input3[site_down][3][n].real());
						projection_spinor_minus2[n] = std::complex<real_t>(input4[site_down][1][n].real() - input4[site_down][3][n].imag(), input4[site_down][1][n].imag() + input4[site_down][3][n].real());
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = std::complex<real_t>(input3[site_up][1][n].real() + input3[site_up][3][n].imag(), input3[site_up][1][n].imag() - input3[site_up][3][n].real());
						projection_spinor_plus2[n] = std::complex<real_t>(input4[site_up][1][n].real() + input4[site_up][3][n].imag(), input4[site_up][1][n].imag() - input4[site_up][3][n].real());
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for(int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						}
						for(int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
						}

						output3[site][1][i] -= kappa*(tmp1+tmm1);
						output4[site][1][i] -= kappa*(tmp2+tmm2);
						output3[site][3][i] += kappa*std::complex<real_t>(tmp1.imag() - tmm1.imag(), tmm1.real() - tmp1.real());
						output4[site][3][i] += kappa*std::complex<real_t>(tmp2.imag() - tmm2.imag(), tmm2.real() - tmp2.real());
					}
				}
			}



			{
				const size_t site_down = Vector::sdn(site,3);
				const size_t site_up = Vector::sup(site,3);
				{

					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input1[site_down][0][n] + (input1[site_down][2][n]);
						projection_spinor_minus2[n] = input2[site_down][0][n] + (input2[site_down][2][n]);
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input1[site_up][0][n] - (input1[site_up][2][n]);
						projection_spinor_plus2[n] = input2[site_up][0][n] - (input2[site_up][2][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						}

						output1[site][0][i] -= kappa*(tmp1+tmm1);
						output2[site][0][i] -= kappa*(tmp2+tmm2);
						output1[site][2][i] += kappa*(tmp1-tmm1);
						output2[site][2][i] += kappa*(tmp2-tmm2);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][3][n]);
						projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][3][n]);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input1[site_up][1][n] - (input1[site_up][3][n]);
						projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][3][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						}

						output1[site][1][i] -= kappa*(tmp1+tmm1);
						output2[site][1][i] -= kappa*(tmp2+tmm2);
						output1[site][3][i] += kappa*(tmp1-tmm1);
						output2[site][3][i] += kappa*(tmp2-tmm2);
					}
				}
				{
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input3[site_down][0][n] + (input3[site_down][2][n]);
						projection_spinor_minus2[n] = input4[site_down][0][n] + (input4[site_down][2][n]);
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input3[site_up][0][n] - (input3[site_up][2][n]);
						projection_spinor_plus2[n] = input4[site_up][0][n] - (input4[site_up][2][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						}

						output3[site][0][i] -= kappa*(tmp1+tmm1);
						output4[site][0][i] -= kappa*(tmp2+tmm2);
						output3[site][2][i] += kappa*(tmp1-tmm1);
						output4[site][2][i] += kappa*(tmp2-tmm2);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_minus1[n] = input3[site_down][1][n] + (input3[site_down][3][n]);
						projection_spinor_minus2[n] = input4[site_down][1][n] + (input4[site_down][3][n]);
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						projection_spinor_plus1[n] = input3[site_up][1][n] - (input3[site_up][3][n]);
						projection_spinor_plus2[n] = input4[site_up][1][n] - (input4[site_up][3][n]);
					}
					for (int i = 0; i < diracVectorLength; ++i) {
						tmp1 = 0;
						tmp2 = 0;
						tmm1 = 0;
						tmm2 = 0;
						for (int n = 0; n < diracVectorLength; ++n) {
							tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
							tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						}
						for (int n = 0; n < diracVectorLength; ++n) {
							tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
							tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
						}

						output3[site][1][i] -= kappa*(tmp1+tmm1);
						output4[site][1][i] -= kappa*(tmp2+tmm2);
						output3[site][3][i] += kappa*(tmp1-tmm1);
						output4[site][3][i] += kappa*(tmp2-tmm2);
					}
				}
			}
		}

	}
	output1.updateHalo();//TODO is needed?
	output2.updateHalo();//TODO is needed?
	output3.updateHalo();//TODO is needed?
	output4.updateHalo();//TODO is needed?
}

void DiracWilsonOperator::multiply(reduced_dirac_vector_t& output1, reduced_dirac_vector_t& output2, const reduced_dirac_vector_t& input1, const reduced_dirac_vector_t& input2) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;
	const reduced_fermion_lattice_t& linkconf = (lattice);

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
		std::complex<real_t> projection_spinor_minus1[diracVectorLength], projection_spinor_plus1[diracVectorLength], tmm1, tmp1;
		std::complex<real_t> projection_spinor_minus2[diracVectorLength], projection_spinor_plus2[diracVectorLength], tmm2, tmp2;
		{
			{
				const size_t site_down = Vector::sdn(site,0);
				const size_t site_up = Vector::sup(site,0);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real()+input1[site_down][3][n].imag(),input1[site_down][0][n].imag()-input1[site_down][3][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real()+input2[site_down][3][n].imag(),input2[site_down][0][n].imag()-input2[site_down][3][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real()-input1[site_up][3][n].imag(),input1[site_up][0][n].imag()+input1[site_up][3][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real()-input2[site_up][3][n].imag(),input2[site_up][0][n].imag()+input2[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp1 = 0;
					tmm2 = 0;
					tmm2 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
					}
					output1[site][0][i] = input1[site][0][i] - kappa*(tmp1+tmm2);
					output2[site][0][i] = input2[site][0][i] - kappa*(tmp2+tmm2);
					output1[site][3][i] = -input1[site][3][i] + kappa*std::complex<real_t>(tmm1.imag() - tmp1.imag(),tmp1.real() - tmm1.real());
					output2[site][3][i] = -input2[site][3][i] + kappa*std::complex<real_t>(tmm2.imag() - tmp2.imag(),tmp2.real() - tmm2.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real()+input1[site_down][2][n].imag(), input1[site_down][1][n].imag()-input1[site_down][2][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real()+input2[site_down][2][n].imag(), input2[site_down][1][n].imag()-input2[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real()-input1[site_up][2][n].imag(),input1[site_up][1][n].imag()+input1[site_up][2][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real()-input2[site_up][2][n].imag(),input2[site_up][1][n].imag()+input2[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,0)][0](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][0](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][0](i,n);
					}

					output1[site][1][i] = input1[site][1][i]- kappa*(tmp1+tmm1);
					output2[site][1][i] = input2[site][1][i]- kappa*(tmp2+tmm2);
					output1[site][2][i] = -input1[site][2][i]+ kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real()-tmm1.real());
					output2[site][2][i] = -input2[site][2][i]+ kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real()-tmm2.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,1);
				const size_t site_up = Vector::sup(site,1);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][0][n] - (input1[site_down][3][n]);
					projection_spinor_minus2[n] = input2[site_down][0][n] - (input2[site_down][3][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][0][n] + (input1[site_up][3][n]);
					projection_spinor_plus2[n] = input2[site_up][0][n] + (input2[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output1[site][3][i] += kappa*(tmm1-tmp1);
					output2[site][3][i] += kappa*(tmm2-tmp2);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][2][n]);
					projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][2][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input2[site_up][1][n] - (input1[site_up][2][n]);
					projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm2 = 0;
					tmm1 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,1)][1](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][1](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][1](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output1[site][2][i] += kappa*(tmp1-tmm1);
					output2[site][2][i] += kappa*(tmp2-tmm2);
				}
			}
			{
				const size_t site_down = Vector::sdn(site,2);
				const size_t site_up = Vector::sup(site,2);
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][0][n].real() + input1[site_down][2][n].imag(), input1[site_down][0][n].imag() - input1[site_down][2][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][0][n].real() + input2[site_down][2][n].imag(), input2[site_down][0][n].imag() - input2[site_down][2][n].real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][0][n].real() - input1[site_up][2][n].imag(), input1[site_up][0][n].imag() + input1[site_up][2][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][0][n].real() - input2[site_up][2][n].imag(), input2[site_up][0][n].imag() + input2[site_up][2][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output1[site][2][i] += kappa*std::complex<real_t>(tmm1.imag()-tmp1.imag(),tmp1.real() - tmm1.real());
					output2[site][2][i] += kappa*std::complex<real_t>(tmm2.imag()-tmp2.imag(),tmp2.real() - tmm2.real());
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = std::complex<real_t>(input1[site_down][1][n].real() - input1[site_down][3][n].imag(), input1[site_down][1][n].imag() + input1[site_down][3][n].real());
					projection_spinor_minus2[n] = std::complex<real_t>(input2[site_down][1][n].real() - input2[site_down][3][n].imag(), input2[site_down][1][n].imag() + input2[site_down][3][n].real());
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = std::complex<real_t>(input1[site_up][1][n].real() + input1[site_up][3][n].imag(), input1[site_up][1][n].imag() - input1[site_up][3][n].real());
					projection_spinor_plus2[n] = std::complex<real_t>(input2[site_up][1][n].real() + input2[site_up][3][n].imag(), input2[site_up][1][n].imag() - input2[site_up][3][n].real());
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for(int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,2)][2](n,i));
					}
					for(int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][2](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][2](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output1[site][3][i] += kappa*std::complex<real_t>(tmp1.imag() - tmm1.imag(), tmm1.real() - tmp1.real());
					output2[site][3][i] += kappa*std::complex<real_t>(tmp2.imag() - tmm2.imag(), tmm2.real() - tmp2.real());
				}
			}
			{
				const size_t site_down = Vector::sdn(site,3);
				const size_t site_up = Vector::sup(site,3);
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][0][n] + (input1[site_down][2][n]);
					projection_spinor_minus2[n] = input2[site_down][0][n] + (input2[site_down][2][n]);
				}
				for(int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][0][n] - (input1[site_up][2][n]);
					projection_spinor_plus2[n] = input2[site_up][0][n] - (input2[site_up][2][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
					}

					output1[site][0][i] -= kappa*(tmp1+tmm1);
					output2[site][0][i] -= kappa*(tmp2+tmm2);
					output1[site][2][i] += kappa*(tmp1-tmm1);
					output2[site][2][i] += kappa*(tmp2-tmm2);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_minus1[n] = input1[site_down][1][n] + (input1[site_down][3][n]);
					projection_spinor_minus2[n] = input2[site_down][1][n] + (input2[site_down][3][n]);
				}
				for (int n = 0; n < diracVectorLength; ++n) {
					projection_spinor_plus1[n] = input1[site_up][1][n] - (input1[site_up][3][n]);
					projection_spinor_plus2[n] = input2[site_up][1][n] - (input2[site_up][3][n]);
				}
				for (int i = 0; i < diracVectorLength; ++i) {
					tmp1 = 0;
					tmp2 = 0;
					tmm1 = 0;
					tmm2 = 0;
					for (int n = 0; n < diracVectorLength; ++n) {
						tmp1 += projection_spinor_minus1[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
						tmp2 += projection_spinor_minus2[n] * conj(linkconf[Lattice::sdn(site,3)][3](n,i));
					}
					for (int n = 0; n < diracVectorLength; ++n) {
						tmm1 += projection_spinor_plus1[n] * linkconf[site][3](i,n);
						tmm2 += projection_spinor_plus2[n] * linkconf[site][3](i,n);
					}

					output1[site][1][i] -= kappa*(tmp1+tmm1);
					output2[site][1][i] -= kappa*(tmp2+tmm2);
					output1[site][3][i] += kappa*(tmp1-tmm1);
					output2[site][3][i] += kappa*(tmp2-tmm2);
				}
			}
		}
	}
	output1.updateHalo();//TODO is needed?
	output2.updateHalo();//TODO is needed?
}

void DiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;
	const reduced_fermion_lattice_t& linkconf = (lattice);

	/*typedef std::complex<double> SiteElem[4][diracVectorLength];

	const std::complex<double>* tmp = reinterpret_cast<const std::complex<double>*>(inputt.getRaw());
	const SiteElem* input = ((const std::complex<double> (*)[4][3])tmp);

	std::complex<double>* tmp2 = reinterpret_cast<std::complex<double>*>(outputt.getRaw());
	SiteElem* output = ((std::complex<double> (*)[4][3])tmp2);*/
	
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

/*void DiracWilsonOperator::multiply(dirac_vector_t& output, const dirac_vector_t& input) {
	typedef adjoint_lattice_t LT;
	typedef adjoint_lattice_t Vector;
#ifdef MULTITHREADING
	#pragma omp parallel for //shared(output, input, kappa, lattice) default(none)
#endif
	const dirac_vector_t& In = input;
	dirac_vector_t& Out = output;
	const adjoint_lattice_t& linkconf = (*lattice);
	real_t Kappa = this->kappa;
	std::complex<real_t> tmp1[3], tmp2[3], tmm, tmp;

	double* llRe = (double*) memalign(32,sizeof(double)*4);
	double* llIm = (double*) memalign(32,sizeof(double)*4);
	//double* data = (double*) memalign(32,sizeof(double)*(12));
	for (unsigned int site = 0; site < lattice->localsize; ++site) {
		output[site][0] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0] + -htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][0] + ((*lattice)[site][1])*input[LT::sup(site, 1)][3] + ((*lattice)[site][2])*input[LT::sup(site, 2)][0] + (I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][0] + -((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + input[site][0];
		output[site][1] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1] + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1] + (I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][1] + -((*lattice)[site][1])*input[LT::sup(site, 1)][2] + ((*lattice)[site][2])*input[LT::sup(site, 2)][1] + (-I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][1] + -((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + input[site][1];
		output[site][2] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2]) + (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + -(kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][1]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][2]) + (-I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][0])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][0]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + -input[site][2];
		output[site][3] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + -(kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3]) + (-I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][0]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][3]) + (I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][1])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][1]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + -input[site][3];

		//First we start the hopping parameter terms
		GaugeVector tmp[4][4];
		//First we put U(x,mu)*input(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_up = LT::sup(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				tmp[mu][nu] = (*lattice)[site][mu]*input[site_up][nu];
			}
		}
		//float* llRe = (float*) memalign(16,sizeof(float)*4);
		//float* llIm = (float*) memalign(16,sizeof(float)*4);
		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_up = LT::sup(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				//tmp[mu][nu] = (*lattice)[site][mu]*input[site_up][nu];
				//First we load the real and the imaginary part of the input
				//Then we load the transpose of the matrix
				//FermionicGroup g = (*lattice)[site][mu];
				//FermionicGroup* t = (FermionicGroup*) memalign(32,sizeof(FermionicGroup)*4);
				//t[0] = (*lattice)[site][mu];
				//const double* data = (*lattice)[site][mu].data();
				//for (int i = 0; i < 9; ++i) data[i] = (*lattice)[site][mu].data()[i];

				__m256d row1 = _mm256_load_pd(mem+(12*4*site+12*mu+(4*0)));
				__m256d row2 = _mm256_load_pd(mem+(12*4*site+12*mu+(4*1)));
				__m256d row3 = _mm256_load_pd(mem+(12*4*site+12*mu+(4*2)));

				//Then we perform the multiplication
				__m256d res1Re = _mm256_mul_pd(row1,_mm256_set1_pd(real(input[site_up][nu][0])));
				__m256d res2Re = _mm256_mul_pd(row2,_mm256_set1_pd(real(input[site_up][nu][1])));
				__m256d res3Re = _mm256_mul_pd(row3,_mm256_set1_pd(real(input[site_up][nu][2])));
				__m256d res1Im = _mm256_mul_pd(row1,_mm256_set1_pd(imag(input[site_up][nu][0])));
				__m256d res2Im = _mm256_mul_pd(row2,_mm256_set1_pd(imag(input[site_up][nu][1])));
				__m256d res3Im = _mm256_mul_pd(row3,_mm256_set1_pd(imag(input[site_up][nu][2])));

				//Now we get the final result
				__m256d resfRe = _mm256_add_pd(_mm256_add_pd(res1Re,res2Re),res3Re);
				__m256d resfIm = _mm256_add_pd(_mm256_add_pd(res1Im,res2Im),res3Im);

				//We get the result
				_mm256_store_pd(llRe,resfRe);
				_mm256_store_pd(llIm,resfIm);
				//Everything is stored in the reverse order!
				tmp[mu][nu][0] = std::complex<real_t>(llRe[3],llIm[3]);
				tmp[mu][nu][1] = std::complex<real_t>(llRe[2],llIm[2]);
				tmp[mu][nu][3] = std::complex<real_t>(llRe[1],llIm[1]);


				//Then we load the transpose of the matrix
				__m128 row1 = _mm_set_ps((*lattice)[site][mu].at(0,0),(*lattice)[site][mu].at(1,0),(*lattice)[site][mu].at(2,0), 0);
				__m128 row2 = _mm_set_ps((*lattice)[site][mu].at(0,1),(*lattice)[site][mu].at(1,1),(*lattice)[site][mu].at(2,1), 0);
				__m128 row3 = _mm_set_ps((*lattice)[site][mu].at(0,2),(*lattice)[site][mu].at(1,2),(*lattice)[site][mu].at(2,2), 0);

				//Then we perform the multiplication
				__m128 res1Re = _mm_mul_ps(row1,_mm_set_ps1(real(input[site_up][nu][0])));
				__m128 res2Re = _mm_mul_ps(row2,_mm_set_ps1(real(input[site_up][nu][1])));
				__m128 res3Re = _mm_mul_ps(row3,_mm_set_ps1(real(input[site_up][nu][2])));
				__m128 res1Im = _mm_mul_ps(row1,_mm_set_ps1(imag(input[site_up][nu][0])));
				__m128 res2Im = _mm_mul_ps(row2,_mm_set_ps1(imag(input[site_up][nu][1])));
				__m128 res3Im = _mm_mul_ps(row3,_mm_set_ps1(imag(input[site_up][nu][2])));

				//Now we get the final result
				//__m128 resfRe = _mm_add_ps(_mm_add_ps(res1Re,res2Re),res3Re);
				//__m128 resfIm = _mm_add_ps(_mm_add_ps(res1Im,res2Im),res3Im);

				//We get the result
				_mm_store_ps(llRe,_mm_add_ps(_mm_add_ps(res1Re,res2Re),res3Re));
				_mm_store_ps(llIm,_mm_add_ps(_mm_add_ps(res1Im,res2Im),res3Im));
				//Everything is stored in the reverse order!
				tmp[mu][nu][0] = std::complex<real_t>(llRe[3],llIm[3]);
				tmp[mu][nu][1] = std::complex<real_t>(llRe[2],llIm[2]);
				tmp[mu][nu][2] = std::complex<real_t>(llRe[1],llIm[1]);
			}
		}
		//free(llRe);
		//free(llIm);
		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4];
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			hopping[0][alpha] = (tmp[0][0][alpha]+I*tmp[0][3][alpha]+tmp[1][0][alpha]+tmp[1][3][alpha]+tmp[2][0][alpha]+I*tmp[2][2][alpha]+tmp[3][0][alpha]-tmp[3][2][alpha]);
			hopping[1][alpha] = (tmp[0][1][alpha]+I*tmp[0][2][alpha]+tmp[1][1][alpha]-tmp[1][2][alpha]+tmp[2][1][alpha]-I*tmp[2][3][alpha]+tmp[3][1][alpha]-tmp[3][3][alpha]);
			hopping[2][alpha] = (I*tmp[0][1][alpha]-tmp[0][2][alpha]+tmp[1][1][alpha]-tmp[1][2][alpha]+I*tmp[2][0][alpha]-tmp[2][2][alpha]+tmp[3][0][alpha]-tmp[3][2][alpha]);
			hopping[3][alpha] = (I*tmp[0][0][alpha]-tmp[0][3][alpha]-tmp[1][0][alpha]-tmp[1][3][alpha]-I*tmp[2][1][alpha]-tmp[2][3][alpha]+tmp[3][1][alpha]-tmp[3][3][alpha]);
		}

		//Then we put U(x-mu,mu)*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			int site_down = LT::sdn(site,mu);
			for (unsigned int nu = 0; nu < 4; ++nu) {
				//tmp[mu][nu] = htrans((*lattice)[site_down][mu])*input[site_down][nu];

				__m256d row1 = _mm256_load_pd(mem+(12*4*site_down+12*mu+(4*0)));
				__m256d row2 = _mm256_load_pd(mem+(12*4*site_down+12*mu+(4*1)));
				__m256d row3 = _mm256_load_pd(mem+(12*4*site_down+12*mu+(4*2)));

				//Then we perform the multiplication
				__m256d res1Re = _mm256_mul_pd(row1,_mm256_set1_pd(real(input[site_down][nu][0])));
				__m256d res2Re = _mm256_mul_pd(row2,_mm256_set1_pd(real(input[site_down][nu][1])));
				__m256d res3Re = _mm256_mul_pd(row3,_mm256_set1_pd(real(input[site_down][nu][2])));
				__m256d res1Im = _mm256_mul_pd(row1,_mm256_set1_pd(imag(input[site_down][nu][0])));
				__m256d res2Im = _mm256_mul_pd(row2,_mm256_set1_pd(imag(input[site_down][nu][1])));
				__m256d res3Im = _mm256_mul_pd(row3,_mm256_set1_pd(imag(input[site_down][nu][2])));

				//Now we get the final result
				__m256d resfRe = _mm256_add_pd(_mm256_add_pd(res1Re,res2Re),res3Re);
				__m256d resfIm = _mm256_add_pd(_mm256_add_pd(res1Im,res2Im),res3Im);

				//We get the result
				_mm256_store_pd(llRe,resfRe);
				_mm256_store_pd(llIm,resfIm);
				//Everything is stored in the reverse order!
				tmp[mu][nu][0] = std::complex<real_t>(llRe[3],llIm[3]);
				tmp[mu][nu][1] = std::complex<real_t>(llRe[2],llIm[2]);
				tmp[mu][nu][3] = std::complex<real_t>(llRe[1],llIm[1]);
			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			hopping[0][alpha] += (tmp[0][0][alpha]-I*tmp[0][3][alpha]+tmp[1][0][alpha]-tmp[1][3][alpha]+tmp[2][0][alpha]-I*tmp[2][2][alpha]+tmp[3][0][alpha]+tmp[3][2][alpha]);
			hopping[1][alpha] += (tmp[0][1][alpha]-I*tmp[0][2][alpha]+tmp[1][1][alpha]+tmp[1][2][alpha]+tmp[2][1][alpha]+I*tmp[2][3][alpha]+tmp[3][1][alpha]+tmp[3][3][alpha]);
			hopping[2][alpha] += (-I*tmp[0][1][alpha]-tmp[0][2][alpha]-tmp[1][1][alpha]-tmp[1][2][alpha]-I*tmp[2][0][alpha]-tmp[2][2][alpha]-tmp[3][0][alpha]-tmp[3][2][alpha]);
			hopping[3][alpha] += (-I*tmp[0][0][alpha]-tmp[0][3][alpha]+tmp[1][0][alpha]-tmp[1][3][alpha]+I*tmp[2][1][alpha]-tmp[2][3][alpha]-tmp[3][1][alpha]-tmp[3][3][alpha]);
		}
		//The final result is gamma5*input - kappa*hopping
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			output[site][0][alpha] = input[site][0][alpha] - kappa*hopping[0][alpha];
			output[site][1][alpha] = input[site][1][alpha] - kappa*hopping[1][alpha];
			output[site][2][alpha] = -(input[site][2][alpha] + kappa*hopping[2][alpha]);
			output[site][3][alpha] = -(input[site][3][alpha] + kappa*hopping[3][alpha]);
		}

		if (site == 5) {
			std::cout << "Nuovo: " << std::endl;
			std::cout << output[site][0] << std::endl;
			std::cout << output[site][1] << std::endl;
			std::cout << output[site][2] << std::endl;
			std::cout << output[site][3] << std::endl;
			std::cout << "Vecchio: " << std::endl;
			output[site][0] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0] + -htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][0] + ((*lattice)[site][1])*input[LT::sup(site, 1)][3] + ((*lattice)[site][2])*input[LT::sup(site, 2)][0] + (I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][0] + -((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + input[site][0];
			output[site][1] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1] + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1] + (I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][1] + -((*lattice)[site][1])*input[LT::sup(site, 1)][2] + ((*lattice)[site][2])*input[LT::sup(site, 2)][1] + (-I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][1] + -((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + input[site][1];
			output[site][2] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2]) + (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + -(kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][1]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][2]) + (-I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][0])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][0]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + -input[site][2];
			output[site][3] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + -(kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3]) + (-I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][0]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][3]) + (I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][1])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][1]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + -input[site][3];
			std::cout << output[site][0] << std::endl;
			std::cout << output[site][1] << std::endl;
			std::cout << output[site][2] << std::endl;
			std::cout << output[site][3] << std::endl;
			exit(0);
		}

		/*{
			{
				const size_t sitd = Vector::sdn(site,0);
				const size_t situ = Vector::sup(site,0);
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][0][nn] +negatI (In[sitd][3][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][0][nn] +positI (In[situ][3][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][0].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][0].at(ii,nn);
					}
					if(0==0){
						Out[site][0][ii] = In[site][0][ii]- Kappa*(tmp+tmm);
						Out[site][3][ii] = -In[site][3][ii]+ Kappa*(positI(tmp) +negatI(tmm));
					} else{
						Out[site][0][ii] -= Kappa*(tmp+tmm);
						Out[site][3][ii] += Kappa*(positI(tmp) +negatI(tmm));
					}
				};
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][1][nn] +negatI (In[sitd][2][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][1][nn] +positI (In[situ][2][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][0].at(nn,ii);
					} for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][0].at(ii,nn);
					}
					if(0==0){
						Out[site][1][ii] = In[site][1][ii]- Kappa*(tmp+tmm);
						Out[site][2][ii] = -In[site][2][ii]+ Kappa*(positI(tmp) +negatI(tmm));
					} else{
						Out[site][1][ii] -= Kappa*(tmp+tmm);
						Out[site][2][ii] += Kappa*(positI(tmp) +negatI(tmm));
					}
				};
			}
			{
				const size_t sitd = Vector::sdn(site,1);
				const size_t situ = Vector::sup(site,1);
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][0][nn] - (In[sitd][3][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][0][nn] + (In[situ][3][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][1].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][1].at(ii,nn);
					}
					if(1==0){
						Out[site][0][ii] = In[site][0][ii]- Kappa*(tmp+tmm);
						Out[site][3][ii] = -In[site][3][ii]+ Kappa*(-(tmp) +(tmm));
					} else{
						Out[site][0][ii] -= Kappa*(tmp+tmm);
						Out[site][3][ii] += Kappa*(-(tmp) +(tmm));
					}
				};
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][1][nn] + (In[sitd][2][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][1][nn] - (In[situ][2][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][1].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][1].at(ii,nn);
					}
					if(1==0){
						Out[site][1][ii] = In[site][1][ii]- Kappa*(tmp+tmm);
						Out[site][2][ii] = -In[site][2][ii]+ Kappa*((tmp) -(tmm));
					} else{
						Out[site][1][ii] -= Kappa*(tmp+tmm);
						Out[site][2][ii] += Kappa*((tmp) -(tmm));
					}
				};
			}
			{
				const size_t sitd = Vector::sdn(site,2);
				const size_t situ = Vector::sup(site,2);
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][0][nn] +negatI (In[sitd][2][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][0][nn] +positI (In[situ][2][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][2].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][2].at(ii,nn);
					}
					if(2==0){
						Out[site][0][ii] = In[site][0][ii]- Kappa*(tmp+tmm);
						Out[site][2][ii] = -In[site][2][ii]+ Kappa*(positI(tmp) +negatI(tmm));
					} else{
						Out[site][0][ii] -= Kappa*(tmp+tmm);
						Out[site][2][ii] += Kappa*(positI(tmp) +negatI(tmm));
					}
				};
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][1][nn] +positI (In[sitd][3][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][1][nn] +negatI (In[situ][3][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][2].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][2].at(ii,nn);
					}
					if(2==0){
						Out[site][1][ii] = In[site][1][ii]- Kappa*(tmp+tmm);
						Out[site][3][ii] = -In[site][3][ii]+ Kappa*(negatI(tmp) +positI(tmm));
					} else{
						Out[site][1][ii] -= Kappa*(tmp+tmm);
						Out[site][3][ii] += Kappa*(negatI(tmp) +positI(tmm));
					}
				};
			}
			{
				const size_t sitd = Vector::sdn(site,3);
				const size_t situ = Vector::sup(site,3);
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][0][nn] + (In[sitd][2][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][0][nn] - (In[situ][2][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][3].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][3].at(ii,nn);
					}
					if(3==0){
						Out[site][0][ii] = In[site][0][ii]- Kappa*(tmp+tmm);
						Out[site][2][ii] = -In[site][2][ii]+ Kappa*((tmp) -(tmm));
					} else{
						Out[site][0][ii] -= Kappa*(tmp+tmm);
						Out[site][2][ii] += Kappa*((tmp) -(tmm));
					}
				};
				for(size_t nn(3); nn--;){
					tmp1[nn] = In[sitd][1][nn] + (In[sitd][3][nn]);
				}
				for(size_t nn(3); nn--;){
					tmp2[nn] = In[situ][1][nn] - (In[situ][3][nn]);
				}
				for(size_t ii(3); ii--;){
					tmp = 0;
					tmm = 0;
					for(size_t nn(3); nn--;){
						tmp += tmp1[nn] * (*lattice)[sitd][3].at(nn,ii);
					}
					for(size_t nn(3); nn--;){
						tmm += tmp2[nn] * (*lattice)[site][3].at(ii,nn);
					}
					if(3==0){
						Out[site][1][ii] = In[site][1][ii]- Kappa*(tmp+tmm);
						Out[site][3][ii] = -In[site][3][ii]+ Kappa*((tmp) -(tmm));
					} else{
						Out[site][1][ii] -= Kappa*(tmp+tmm);
						Out[site][3][ii] += Kappa*((tmp) -(tmm));
					}
				};
			}
		};

	}

	free(llRe);
	free(llIm);
	output.updateHalo();//TODO is needed?
}*/

/*void DiracWilsonOperator::multiply(dirac_vector_t& output, const dirac_vector_t& input) {
	typedef adjoint_lattice_t LT;
	for (unsigned int site = 0; site < lattice->localsize; site += 4) {
		//First we start the hopping parameter terms
		__m128 tmpReal[4][4][3];
		__m128 tmpImag[4][4][3];

		//First we put U(x,mu)*input(x+mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				for (unsigned int i = 0; i < 3; ++i) {
					__m128 colmat = _mm_set_ps((*lattice)[site][mu].at(i,0),(*lattice)[site+1][mu].at(i,0),(*lattice)[site+2][mu].at(i,0),(*lattice)[site+3][mu].at(i,0));
					__m128 realvect = _mm_set_ps(real(input[LT::sup(site,mu)][nu][0]),real(input[LT::sup(site+1,mu)][nu][0]),real(input[LT::sup(site+2,mu)][nu][0]),real(input[LT::sup(site+3,mu)][nu][0]));
					__m128 imagvect = _mm_set_ps(imag(input[LT::sup(site,mu)][nu][0]),imag(input[LT::sup(site+1,mu)][nu][0]),imag(input[LT::sup(site+2,mu)][nu][0]),imag(input[LT::sup(site+3,mu)][nu][0]));
					tmpReal[mu][nu][i] = _mm_mul_ps(colmat,realvect);
					tmpImag[mu][nu][i] = _mm_mul_ps(colmat,imagvect);
					for (unsigned int j = 1; j < 3; ++j) {
						colmat = _mm_set_ps((*lattice)[site][mu].at(i,j),(*lattice)[site+1][mu].at(i,j),(*lattice)[site+2][mu].at(i,j),(*lattice)[site+3][mu].at(i,j));
						realvect = _mm_set_ps(real(input[LT::sup(site,mu)][nu][j]),real(input[LT::sup(site+1,mu)][nu][j]),real(input[LT::sup(site+2,mu)][nu][j]),real(input[LT::sup(site+3,mu)][nu][j]));
						imagvect = _mm_set_ps(imag(input[LT::sup(site,mu)][nu][j]),imag(input[LT::sup(site+1,mu)][nu][j]),imag(input[LT::sup(site+2,mu)][nu][j]),imag(input[LT::sup(site+3,mu)][nu][j]));
						tmpReal[mu][nu][i] = _mm_add_ps(tmpReal[mu][nu][i],_mm_mul_ps(colmat,realvect));
						tmpImag[mu][nu][i] = _mm_add_ps(tmpImag[mu][nu][i],_mm_mul_ps(colmat,imagvect));
					}
				}

			}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		__m128 hoppingReal[4][3];
		__m128 hoppingImag[4][3];
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			hoppingReal[0][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpReal[0][0][alpha],tmpReal[1][0][alpha]),tmpReal[1][3][alpha]),tmpReal[2][0][alpha]),tmpReal[3][0][alpha]),tmpImag[0][3][alpha]),tmpImag[2][2][alpha]),tmpReal[3][2][alpha]);
			hoppingReal[1][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[2][3][alpha],tmpReal[0][1][alpha]),tmpReal[1][1][alpha]),tmpReal[2][1][alpha]),tmpReal[3][1][alpha]),tmpImag[0][2][alpha]),tmpReal[1][2][alpha]),tmpReal[3][3][alpha]);
			hoppingReal[2][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpReal[1][1][alpha],tmpReal[3][0][alpha]),tmpImag[0][1][alpha]),tmpImag[2][0][alpha]),tmpReal[0][2][alpha]),tmpReal[1][2][alpha]),tmpReal[2][2][alpha]),tmpReal[3][2][alpha]);
			hoppingReal[3][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpImag[2][1][alpha],tmpReal[3][1][alpha]),tmpImag[0][0][alpha]),tmpReal[0][3][alpha]),tmpReal[1][0][alpha]),tmpReal[1][3][alpha]),tmpReal[2][3][alpha]),tmpReal[3][3][alpha]);

			hoppingImag[0][alpha] = _mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][0][alpha],tmpImag[1][0][alpha]),tmpImag[1][3][alpha]),tmpImag[2][0][alpha]),tmpImag[3][0][alpha]),tmpReal[0][3][alpha]),tmpReal[2][2][alpha]),tmpImag[3][2][alpha]);
			hoppingImag[1][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][1][alpha],tmpImag[1][1][alpha]),tmpImag[2][1][alpha]),tmpImag[3][1][alpha]),tmpReal[0][2][alpha]),tmpImag[1][2][alpha]),tmpImag[3][3][alpha]),tmpReal[2][3][alpha]);
			hoppingImag[2][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[1][1][alpha],tmpImag[3][0][alpha]),tmpReal[0][1][alpha]),tmpReal[2][0][alpha]),tmpImag[0][2][alpha]),tmpImag[1][2][alpha]),tmpImag[2][2][alpha]),tmpImag[3][2][alpha]);
			hoppingImag[3][alpha] = _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpImag[3][1][alpha],tmpReal[0][0][alpha]),tmpImag[0][3][alpha]),tmpImag[1][0][alpha]),tmpImag[1][3][alpha]),tmpImag[2][3][alpha]),tmpImag[3][3][alpha]),tmpReal[2][1][alpha]);
		}

		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				for (unsigned int i = 0; i < 3; ++i) {
					__m128 colmat = _mm_set_ps((*lattice)[LT::sdn(site,mu)][mu].at(0,i),(*lattice)[LT::sdn(site+1,mu)][mu].at(0,i),(*lattice)[LT::sdn(site+2,mu)][mu].at(0,i),(*lattice)[LT::sdn(site+3,mu)][mu].at(0,i));
					__m128 realvect = _mm_set_ps(real(input[LT::sdn(site,mu)][nu][0]),real(input[LT::sdn(site+1,mu)][nu][0]),real(input[LT::sdn(site+2,mu)][nu][0]),real(input[LT::sdn(site+3,mu)][nu][0]));
					__m128 imagvect = _mm_set_ps(imag(input[LT::sdn(site,mu)][nu][0]),imag(input[LT::sdn(site+1,mu)][nu][0]),imag(input[LT::sdn(site+2,mu)][nu][0]),imag(input[LT::sdn(site+3,mu)][nu][0]));
					tmpReal[mu][nu][i] = _mm_mul_ps(colmat,realvect);
					tmpImag[mu][nu][i] = _mm_mul_ps(colmat,imagvect);
					for (unsigned int j = 1; j < 3; ++j) {
						colmat = _mm_set_ps((*lattice)[site][mu].at(j,i),(*lattice)[site+1][mu].at(j,i),(*lattice)[site+2][mu].at(j,i),(*lattice)[site+3][mu].at(j,i));
						realvect = _mm_set_ps(real(input[LT::sdn(site,mu)][nu][j]),real(input[LT::sdn(site+1,mu)][nu][j]),real(input[LT::sdn(site+2,mu)][nu][j]),real(input[LT::sdn(site+3,mu)][nu][j]));
						imagvect = _mm_set_ps(imag(input[LT::sdn(site,mu)][nu][j]),imag(input[LT::sdn(site+1,mu)][nu][j]),imag(input[LT::sdn(site+2,mu)][nu][j]),imag(input[LT::sdn(site+3,mu)][nu][j]));
						tmpReal[mu][nu][i] = _mm_add_ps(tmpReal[mu][nu][i],_mm_mul_ps(colmat,realvect));
						tmpImag[mu][nu][i] = _mm_add_ps(tmpImag[mu][nu][i],_mm_mul_ps(colmat,imagvect));
					}
				}

			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			hoppingReal[0][alpha] = _mm_add_ps(hoppingReal[0][alpha], _mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][3][alpha],tmpImag[2][2][alpha]),tmpReal[0][0][alpha]),tmpReal[1][0][alpha]),tmpReal[2][0][alpha]),tmpReal[3][0][alpha]),tmpReal[3][2][alpha]),tmpReal[1][3][alpha]));
			hoppingReal[1][alpha] = _mm_add_ps(hoppingReal[1][alpha], _mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][2][alpha],tmpReal[0][1][alpha]),tmpReal[1][1][alpha]),tmpReal[1][2][alpha]),tmpReal[2][1][alpha]),tmpReal[3][1][alpha]),tmpReal[3][3][alpha]),tmpImag[2][3][alpha]));
			hoppingReal[2][alpha] = _mm_add_ps(hoppingReal[2][alpha], _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpImag[0][1][alpha],tmpImag[2][0][alpha]),tmpReal[0][2][alpha]),tmpReal[1][1][alpha]),tmpReal[1][2][alpha]),tmpReal[2][2][alpha]),tmpReal[3][0][alpha]),tmpReal[3][2][alpha]));
			hoppingReal[3][alpha] = _mm_add_ps(hoppingReal[3][alpha], _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpImag[0][1][alpha],tmpImag[2][0][alpha]),tmpReal[0][2][alpha]),tmpReal[1][1][alpha]),tmpReal[1][2][alpha]),tmpReal[2][2][alpha]),tmpReal[3][0][alpha]),tmpReal[3][2][alpha]));


			hoppingImag[0][alpha] = _mm_add_ps(hoppingImag[0][alpha], _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][0][alpha],tmpImag[1][0][alpha]),tmpImag[2][0][alpha]),tmpImag[3][0][alpha]),tmpImag[3][2][alpha]),tmpImag[1][3][alpha]),tmpReal[0][3][alpha]),tmpReal[2][2][alpha]));
			hoppingImag[1][alpha] = _mm_add_ps(hoppingImag[1][alpha], _mm_sub_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][1][alpha],tmpImag[1][1][alpha]),tmpImag[1][2][alpha]),tmpImag[2][1][alpha]),tmpImag[3][1][alpha]),tmpImag[3][3][alpha]),tmpReal[2][3][alpha]),tmpReal[0][2][alpha]));
			//Exception!!!
			hoppingImag[2][alpha] = _mm_sub_ps(hoppingImag[2][alpha],_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(_mm_add_ps(tmpImag[0][2][alpha],tmpImag[1][1][alpha]),tmpImag[1][2][alpha]),tmpImag[2][2][alpha]),tmpImag[3][0][alpha]),tmpImag[3][2][alpha]),tmpReal[0][1][alpha]),tmpReal[2][0][alpha]));
			hoppingImag[3][alpha] = _mm_add_ps(hoppingImag[3][alpha], _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_sub_ps(_mm_add_ps(tmpImag[1][0][alpha],tmpReal[2][1][alpha]),tmpImag[0][3][alpha]),tmpImag[1][3][alpha]),tmpImag[2][3][alpha]),tmpImag[3][1][alpha]),tmpImag[3][3][alpha]),tmpReal[0][0][alpha]));
		}

		//Now we multiply by kappa
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int alpha = 0; alpha < 3; ++alpha) {
				hoppingReal[mu][alpha] = _mm_mul_ps(hoppingReal[mu][alpha],_mm_set_ps1(kappa));
				hoppingImag[mu][alpha] = _mm_mul_ps(hoppingImag[mu][alpha],_mm_set_ps1(kappa));
			}
		}

		float* memR = (float*) memalign(16,sizeof(float)*4);
		float* memI = (float*) memalign(16,sizeof(float)*4);

		//The final result is gamma5*input - kappa*hopping
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			for (unsigned int delta = 0; delta < 4; ++delta) {
				output[site+delta][0][alpha] = input[site+delta][0][alpha];
				output[site+delta][1][alpha] = input[site+delta][1][alpha];
				output[site+delta][2][alpha] = -input[site+delta][2][alpha];
				output[site+delta][3][alpha] = -input[site+delta][3][alpha];
			}
		}
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				_mm_storer_ps(memR,hoppingReal[mu][alpha]);
				_mm_storer_ps(memI,hoppingImag[mu][alpha]);
				for (unsigned int delta = 0; delta < 4; ++delta) {
					output[site+delta][mu][alpha] -= std::complex<real_t>(memR[delta],memI[delta]);
				}
			}
		}

		free(memR);
		free(memI);
	}

	{
		int site = 5;
		std::cout << "Nuovo: " << std::endl;
		std::cout << output[site][0] << std::endl;
		std::cout << output[site][1] << std::endl;
		std::cout << output[site][2] << std::endl;
		std::cout << output[site][3] << std::endl;
		std::cout << "Vecchio: " << std::endl;
		output[site][0] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0] + -htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][0] + ((*lattice)[site][1])*input[LT::sup(site, 1)][3] + ((*lattice)[site][2])*input[LT::sup(site, 2)][0] + (I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][0] + -((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + input[site][0];
		output[site][1] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1] + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1] + (I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][1] + -((*lattice)[site][1])*input[LT::sup(site, 1)][2] + ((*lattice)[site][2])*input[LT::sup(site, 2)][1] + (-I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][1] + -((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + input[site][1];
		output[site][2] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2]) + (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + -(kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][1]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][2]) + (-I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][0])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][0]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + -input[site][2];
		output[site][3] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + -(kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3]) + (-I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][0]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][3]) + (I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][1])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][1]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + -input[site][3];
		std::cout << output[site][0] << std::endl;
		std::cout << output[site][1] << std::endl;
		std::cout << output[site][2] << std::endl;
		std::cout << output[site][3] << std::endl;
		exit(0);
	}

	output.updateHalo();//TODO is needed?
}*/
/*#define NN 4
void DiracWilsonOperator::multiply(dirac_vector_t& output, const dirac_vector_t& input) {
	typedef adjoint_lattice_t LT;
	typedef adjoint_lattice_t Vector;
#ifdef MULTITHREADING
	#pragma omp parallel for //shared(output, input, kappa, lattice) default(none)
#endif
	for (unsigned int site = 0; site < lattice->localsize; site += NN) {
		//First we start the hopping parameter terms
		std::complex<double> tmp[4][4][3][NN];
		//First we put U(x,mu)*input(x+mu)
		const double* matrix[NN];
		const std::complex<double>* vect[NN];
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {

				for (unsigned int i = 0; i < NN; ++i) {
					matrix[i] = (*lattice)[site+i][mu].data();
					vect[i] = input[LT::sup(site+i,mu)][nu].data();
				}
				for (unsigned int i = 0; i < 3; ++i) {
					for (unsigned int k = 0; k < NN; ++k) {
						tmp[mu][nu][i][k] = matrix[k][3*i]*vect[k][0];
					}
					for (unsigned int j = 1; j < 3; ++j) {
						for (unsigned int k = 0; k < NN; ++k) {
							tmp[mu][nu][i][k] += matrix[k][3*i+j]*vect[k][j];
						}
					}
				}
			}
		}

		//We store the result in a cache multiplied by gamma5(id-gamma[mu]) in dirac space
		GaugeVector hopping[4][4];
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			for (unsigned int delta = 0; delta < NN; ++delta) {
				hopping[0][delta][alpha] = (tmp[0][0][alpha][delta]+I*tmp[0][3][alpha][delta]+tmp[1][0][alpha][delta]+tmp[1][3][alpha][delta]+tmp[2][0][alpha][delta]+I*tmp[2][2][alpha][delta]+tmp[3][0][alpha][delta]-tmp[3][2][alpha][delta]);
				hopping[1][delta][alpha] = (tmp[0][1][alpha][delta]+I*tmp[0][2][alpha][delta]+tmp[1][1][alpha][delta]-tmp[1][2][alpha][delta]+tmp[2][1][alpha][delta]-I*tmp[2][3][alpha][delta]+tmp[3][1][alpha][delta]-tmp[3][3][alpha][delta]);
				hopping[2][delta][alpha] = (I*tmp[0][1][alpha][delta]-tmp[0][2][alpha][delta]+tmp[1][1][alpha][delta]-tmp[1][2][alpha][delta]+I*tmp[2][0][alpha][delta]-tmp[2][2][alpha][delta]+tmp[3][0][alpha][delta]-tmp[3][2][alpha][delta]);
				hopping[3][delta][alpha] = (I*tmp[0][0][alpha][delta]-tmp[0][3][alpha][delta]-tmp[1][0][alpha][delta]-tmp[1][3][alpha][delta]-I*tmp[2][1][alpha][delta]-tmp[2][3][alpha][delta]+tmp[3][1][alpha][delta]-tmp[3][3][alpha][delta]);
			}
		}

		//Then we put U(x-mu,mu)^\dag*input(x-mu)
		for (unsigned int mu = 0; mu < 4; ++mu) {
			for (unsigned int nu = 0; nu < 4; ++nu) {
				//tmp[mu][nu][0] = htrans((*lattice)[LT::sdn(site,mu)][mu])*input[LT::sdn(site,mu)][nu];
				//tmp[mu][nu][1] = htrans((*lattice)[LT::sdn(site+1,mu)][mu])*input[LT::sdn(site+1,mu)][nu];
				const double* matrix0 = (*lattice)[LT::sdn(site,mu)][mu].data();
				const double* matrix1 = (*lattice)[LT::sdn(site+1,mu)][mu].data();
				const double* matrix2 = (*lattice)[LT::sdn(site+2,mu)][mu].data();
				const double* matrix3 = (*lattice)[LT::sdn(site+3,mu)][mu].data();
				const std::complex<double>* vect0 = input[LT::sdn(site,mu)][nu].data();
				const std::complex<double>* vect1 = input[LT::sdn(site+1,mu)][nu].data();
				const std::complex<double>* vect2 = input[LT::sdn(site+2,mu)][nu].data();
				const std::complex<double>* vect3 = input[LT::sdn(site+3,mu)][nu].data();
				for (unsigned int i = 0; i < 3; ++i) {
					tmp[mu][nu][i][0] = matrix0[i]*vect0[0];
					tmp[mu][nu][i][1] = matrix1[i]*vect1[0];
					tmp[mu][nu][i][2] = matrix2[i]*vect2[0];
					tmp[mu][nu][i][3] = matrix3[i]*vect3[0];
					for (unsigned int j = 1; j < 3; ++j) {
						tmp[mu][nu][i][0] += matrix0[3*j+i]*vect0[j];
						tmp[mu][nu][i][1] += matrix1[3*j+i]*vect1[j];
						tmp[mu][nu][i][2] += matrix2[3*j+i]*vect2[j];
						tmp[mu][nu][i][3] += matrix3[3*j+i]*vect3[j];
					}
				}


				for (unsigned int i = 0; i < NN; ++i) {
					matrix[i] = (*lattice)[LT::sdn(site+i,mu)][mu].data();
					vect[i] = input[LT::sdn(site+i,mu)][nu].data();
				}
				for (unsigned int i = 0; i < 3; ++i) {
					for (unsigned int k = 0; k < NN; ++k) {
						tmp[mu][nu][i][k] = matrix[k][i]*vect[k][0];
					}
					for (unsigned int j = 1; j < 3; ++j) {
						for (unsigned int k = 0; k < NN; ++k) {
							tmp[mu][nu][i][k] += matrix[k][3*j+i]*vect[k][j];
						}
					}
				}
			}
		}

		//We store the result in the same a cache multiplied by gamma5(id+gamma[mu]) in dirac space
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			for (unsigned int delta = 0; delta < NN; ++delta) {
				hopping[0][delta][alpha] += (tmp[0][0][alpha][delta]-I*tmp[0][3][alpha][delta]+tmp[1][0][alpha][delta]-tmp[1][3][alpha][delta]+tmp[2][0][alpha][delta]-I*tmp[2][2][alpha][delta]+tmp[3][0][alpha][delta]+tmp[3][2][alpha][delta]);
				hopping[1][delta][alpha] += (tmp[0][1][alpha][delta]-I*tmp[0][2][alpha][delta]+tmp[1][1][alpha][delta]+tmp[1][2][alpha][delta]+tmp[2][1][alpha][delta]+I*tmp[2][3][alpha][delta]+tmp[3][1][alpha][delta]+tmp[3][3][alpha][delta]);
				hopping[2][delta][alpha] += (-I*tmp[0][1][alpha][delta]-tmp[0][2][alpha][delta]-tmp[1][1][alpha][delta]-tmp[1][2][alpha][delta]-I*tmp[2][0][alpha][delta]-tmp[2][2][alpha][delta]-tmp[3][0][alpha][delta]-tmp[3][2][alpha][delta]);
				hopping[3][delta][alpha] += (-I*tmp[0][0][alpha][delta]-tmp[0][3][alpha][delta]+tmp[1][0][alpha][delta]-tmp[1][3][alpha][delta]+I*tmp[2][1][alpha][delta]-tmp[2][3][alpha][delta]-tmp[3][1][alpha][delta]-tmp[3][3][alpha][delta]);
			}
		}

		//The final result is gamma5*input - kappa*hopping
		for (unsigned int alpha = 0; alpha < 3; ++alpha) {
			for (unsigned int delta = 0; delta < NN; ++delta) {
				output[site+delta][0][alpha] = input[site+delta][0][alpha] - kappa*hopping[0][delta][alpha];
				output[site+delta][1][alpha] = input[site+delta][1][alpha] - kappa*hopping[1][delta][alpha];
				output[site+delta][2][alpha] = -(input[site+delta][2][alpha] + kappa*hopping[2][delta][alpha]);
				output[site+delta][3][alpha] = -(input[site+delta][3][alpha] + kappa*hopping[3][delta][alpha]);
			}
		}

		if (site == 5) {
			std::cout << "Nuovo: " << std::endl;
			std::cout << output[site][0] << std::endl;
			std::cout << output[site][1] << std::endl;
			std::cout << output[site][2] << std::endl;
			std::cout << output[site][3] << std::endl;
			std::cout << "Vecchio: " << std::endl;
			output[site][0] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0] + -htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0] + (-I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][0] + ((*lattice)[site][1])*input[LT::sup(site, 1)][3] + ((*lattice)[site][2])*input[LT::sup(site, 2)][0] + (I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][0] + -((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + input[site][0];
			output[site][1] = -(kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1] + (-I)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1] + htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2] + htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1] + (I)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1] + htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + -(kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1] + (I)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + ((*lattice)[site][1])*input[LT::sup(site, 1)][1] + -((*lattice)[site][1])*input[LT::sup(site, 1)][2] + ((*lattice)[site][2])*input[LT::sup(site, 2)][1] + (-I)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + ((*lattice)[site][3])*input[LT::sup(site, 3)][1] + -((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + input[site][1];
			output[site][2] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][2]) + (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][2]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][2]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][1])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][2]) + -(kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][1]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][2]) + (-I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][0])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][2]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][0]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][2]) + -input[site][2];
			output[site][3] = (I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][0])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 0)][0])*input[LT::sdn(site, 0)][3]) + -(kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][0]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 1)][1])*input[LT::sdn(site, 1)][3]) + (-I)*((kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][1])) + (kappa)*(htrans((*lattice)[LT::sdn(site, 2)][2])*input[LT::sdn(site, 2)][3]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][1]) + (kappa)*(htrans((*lattice)[LT::sdn(site, 3)][3])*input[LT::sdn(site, 3)][3]) + (-I)*((kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][0])) + (kappa)*(((*lattice)[site][0])*input[LT::sup(site, 0)][3]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][0]) + (kappa)*(((*lattice)[site][1])*input[LT::sup(site, 1)][3]) + (I)*((kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][1])) + (kappa)*(((*lattice)[site][2])*input[LT::sup(site, 2)][3]) + -(kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][1]) + (kappa)*(((*lattice)[site][3])*input[LT::sup(site, 3)][3]) + -input[site][3];
			std::cout << output[site][0] << std::endl;
			std::cout << output[site][1] << std::endl;
			std::cout << output[site][2] << std::endl;
			std::cout << output[site][3] << std::endl;
			exit(0);
		}

	}
	output.updateHalo();//TODO is needed?
}*/

 void DiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
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
		 }
		 if (!gamma5) {
			for (int i = 0; i < diracVectorLength; ++i) {
				output[site][2][i] = -output[site][2][i]+static_cast<real_t>(2)*alpha*vector2[site][2][i];
				output[site][3][i] = -output[site][3][i]+static_cast<real_t>(2)*alpha*vector2[site][3][i];
			}
		}
	 }
	 output.updateHalo();//TODO is needed?
}

FermionForce* DiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

} /* namespace Update */
