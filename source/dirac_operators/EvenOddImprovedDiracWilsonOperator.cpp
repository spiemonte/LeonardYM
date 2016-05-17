/*
 * ImprovedDiracWilsonOperator.cpp
 *
 *  Created on: May 4, 2012
 *      Author: spiem_01
 */

#include "EvenOddImprovedDiracWilsonOperator.h"
#include "utils/Gamma.h"

namespace Update {

int fieldIndex[4][4] =  {{0,0,1,2},{0,0,3,4},{0,0,0,5},{0,0,0,0}};

inline real_t conj(const real_t& t) {
	return t;
}

inline std::complex<real_t> multiply_by_I(const std::complex<real_t>& a) {
	return std::complex<real_t>(-imag(a),real(a));
}

EvenOddImprovedDiracWilsonOperator::EvenOddImprovedDiracWilsonOperator() : ImprovedDiracWilsonOperator(), cloverMatrixInverse(0) { }

EvenOddImprovedDiracWilsonOperator::EvenOddImprovedDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, double _kappa, double _csw, bool _gamma5) : ImprovedDiracWilsonOperator(_lattice, _kappa, _csw, _gamma5), cloverMatrixInverse(0) {
	this->calculateInverseEvenEven();
}

EvenOddImprovedDiracWilsonOperator::~EvenOddImprovedDiracWilsonOperator() {
	if (cloverMatrixInverse) delete[] cloverMatrixInverse;
}

void EvenOddImprovedDiracWilsonOperator::multiply(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input) {
	this->multiplyEvenOdd(output, input, EVEN);
	this->multiplyEvenEvenInverse(output);
	this->multiplyEvenOdd(output, output, ODD);
	//this->multiplyOddOdd(output, output, ODD);
	this->multiplyOddOddMinusIdentity(output, input, ODD);

	typedef reduced_dirac_vector_t::Layout Layout;
	
	
#pragma omp parallel for
	for (int site = 0; site< Layout::completesize; ++site) {//Odd part?
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == 0) {
			for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] = input[site][mu];
		}
		else if (gamma5) {
			for (unsigned int mu = 2; mu < 4; ++mu) output[site][mu] = -output[site][mu];
		}
	}
	
}

void EvenOddImprovedDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t & , const reduced_dirac_vector_t &  , const reduced_dirac_vector_t & , const complex& ) {
	if (isOutputProcess()) std::cout << "EvenOddImprovedDiracWilsonOperator::multiplyAdd not implemented" << std::endl;
	exit(5);
}

void EvenOddImprovedDiracWilsonOperator::setLattice(const extended_fermion_lattice_t& _lattice) {
	lattice = _lattice;
	this->updateFieldStrength(_lattice);
	this->calculateInverseEvenEven();
}

//This part of code multiplies R_ee^-1 to the even part of Out and store the result in the even part of Out
void EvenOddImprovedDiracWilsonOperator::multiplyEvenEvenInverse(reduced_dirac_vector_t & output) {
	typedef reduced_dirac_vector_t Vector;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (int site = 0; site< Layout::completesize; ++site) {//Even part?
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == 0) {
			//TODO antialias?
			std::complex<real_t> antialias[4][diracVectorLength];
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				for (unsigned int c = 0; c < diracVectorLength; c++) {
					antialias[alpha][c] = output[site][alpha][c];
				}
			}
			//End antialias
			for (unsigned int alpha = 0; alpha < 4; ++alpha) {
				for (unsigned int c = 0; c < diracVectorLength; c++) {
					std::complex<real_t> tmp = std::complex<real_t>(0.,0.);
					for (unsigned int beta = 0; beta < 4; ++beta) {//sum over beta
						for (unsigned int b = 0; b < diracVectorLength; ++b) {
							tmp += cloverMatrixInverse[site].at(3*alpha+c,3*beta+b) * antialias[beta][b];//TODO disperazione
			    			}
					}
					output[site][alpha][c] = /*complex<Real>(0.,kappa*cloverTerm/2.)**/tmp;//Factor kappa and 2 and inverse ERRORE
				}
			}
		}
	}
}

//This part of code multiplies R_oo to the odd part of In and -D_oe to the even part of Out and store the result in the odd part of Out
void EvenOddImprovedDiracWilsonOperator::multiplyOddOdd(reduced_dirac_vector_t & output, Part part) {
	typedef reduced_dirac_vector_t Vector;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (int site = 0; site< Layout::completesize; ++site) {//Odd part?
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == (int)(part)) {
			//We store the result of the clover term in an intermediate vector
			GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				clover[0][i] = 0;
				clover[1][i] = 0;
				clover[2][i] = 0;
				clover[3][i] = 0;
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += multiply_by_I((-F[site][0].at(i,j)+F[site][5].at(i,j))*output[site][0][j]);
					clover[1][i] += multiply_by_I((+F[site][0].at(i,j)-F[site][5].at(i,j))*output[site][1][j]);
					clover[2][i] += multiply_by_I((+F[site][0].at(i,j)+F[site][5].at(i,j))*output[site][2][j]);
					clover[3][i] += multiply_by_I((-F[site][0].at(i,j)-F[site][5].at(i,j))*output[site][3][j]);
				}
			}

#ifdef ADJOINT
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += std::complex<real_t>(+F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*output[site][1][j];
					clover[1][i] += std::complex<real_t>(-F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)))*output[site][0][j];
					clover[2][i] += std::complex<real_t>(-F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*output[site][3][j];
					clover[3][i] += std::complex<real_t>(+F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)))*output[site][2][j];
				}
			}
#endif
#ifndef ADJOINT
			//We store the result of the clover term in an intermediate vector
			//GaugeVector clover[4];
			for (int i = 0; i < diracVectorLength; ++i) {
				for (int j = 0; j < diracVectorLength; ++j) {
					clover[0][i] += (+F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*output[site][1][j];
					clover[1][i] += (-F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)))*output[site][0][j];
					clover[2][i] += (-F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*output[site][3][j];
					clover[3][i] += (+F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)))*output[site][2][j];
				}
			}
#endif
			
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = output[site][0][n] +(kappa*csw)*clover[0][n];
				output[site][1][n] = output[site][1][n] +(kappa*csw)*clover[1][n];
				output[site][2][n] = output[site][2][n] -(kappa*csw)*clover[2][n];
				output[site][3][n] = output[site][3][n] -(kappa*csw)*clover[3][n];
			}
		}
	}
}

//This part of code multiplies R_oo to the odd part of In and -D_oe to the even part of Out and store the result in the odd part of Out
void EvenOddImprovedDiracWilsonOperator::multiplyOddOddMinusIdentity(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input, Part part) {
	typedef reduced_dirac_vector_t Vector;
	typedef reduced_dirac_vector_t::Layout Layout;

#pragma omp parallel for
	for (int site = 0; site< Layout::completesize; ++site) {//Odd part?
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == (int)(part)) {
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
				output[site][0][n] = input[site][0][n] +(kappa*csw)*clover[0][n] - output[site][0][n];
				output[site][1][n] = input[site][1][n] +(kappa*csw)*clover[1][n] - output[site][1][n];
				output[site][2][n] = input[site][2][n] -(kappa*csw)*clover[2][n]  - output[site][2][n];
				output[site][3][n] = input[site][3][n] -(kappa*csw)*clover[3][n]  - output[site][3][n];
			}
		}
	}
}

//This part of code multiplies D_eo to the odd part of In and store the result in the even part of Out
void EvenOddImprovedDiracWilsonOperator::multiplyEvenOdd(reduced_dirac_vector_t & output, const reduced_dirac_vector_t & input, Part part) {
	typedef reduced_fermion_lattice_t::Layout Layout;
	typedef reduced_dirac_vector_t Vector;

#pragma omp parallel for
	for (int site = 0; site< Layout::localsize; ++site) {//Even part?
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == part) {
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
					GaugeVector tmp = htrans(lattice[Vector::sdn(site,mu)][mu])*projection_spinor[mu][nu];
					tmp_minus[mu][nu] = tmp_plus[mu][nu] - tmp;
					tmp_plus[mu][nu] += tmp;
				}
			}

		
			
				//The final result is - kappa*gamma5*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = std::complex<real_t>( - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])), - kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = std::complex<real_t>( - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])), - kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = std::complex<real_t>( + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])), + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
				output[site][3][n] = std::complex<real_t>( + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])), + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
			}
		}
		else {
			for (unsigned int mu = 0; mu < 4; ++mu) output[site][mu] = input[site][mu];
		}
	}
	output.updateHalo();
}



void EvenOddImprovedDiracWilsonOperator::calculateInverseEvenEven() {
	typedef reduced_fermion_lattice_t::Layout Layout;
	if (cloverMatrixInverse == 0) cloverMatrixInverse = new clover_matrix_t[Layout::completesize];

#pragma omp parallel for
	for (int site = 0; site < Layout::completesize; ++site) {
		/*for (unsigned int alpha = 0; alpha < 4; ++alpha) {
			for (unsigned int beta = 0; beta < 4; ++beta) {
				for (int a = 0; a < diracVectorLength; ++a) {
					for (size_t b = 0; b < diracVectorLength; ++b) {
						if (3*alpha + a == 3*beta +b)  cloverMatrixInverse[site].at(3*alpha + a, 3*beta +b) = std::complex<real_t>(1.,0.);
						else cloverMatrixInverse[site].at(3*alpha + a, 3*beta +b) = std::complex<real_t>(0.,0.);
						for (unsigned int mu = 0; mu < 4; ++mu) {
							for (unsigned int nu = mu + 1; nu < 4; ++nu) {
								//sigma_{\mu\nu} F^{\mu\nu} In
								cloverMatrixInverse[site].at(3*alpha + a, 3*beta +b) += std::complex<real_t>(0.,kappa*csw) * Sigma::sigma(mu, nu, alpha, beta)  * F[site][fieldIndex[mu][nu]].at(a,b);
				    			}
						}
					}
				}
			}
		}*/

		cloverMatrixInverse[site] = clover_matrix_t::Identity();

		for (int i = 0; i < diracVectorLength; ++i) {
			for (int j = 0; j < diracVectorLength; ++j) {
				cloverMatrixInverse[site].at(3*0 + i, 3*0 +j) += (kappa*csw)*std::complex<real_t>(0,1)*(-F[site][0].at(i,j)+F[site][5].at(i,j));
				cloverMatrixInverse[site].at(3*1 + i, 3*1 +j) += (kappa*csw)*std::complex<real_t>(0,1)*(+F[site][0].at(i,j)-F[site][5].at(i,j));
				cloverMatrixInverse[site].at(3*2 + i, 3*2 +j) -= (kappa*csw)*std::complex<real_t>(0,1)*(+F[site][0].at(i,j)+F[site][5].at(i,j));
				cloverMatrixInverse[site].at(3*3 + i, 3*3 +j) -= (kappa*csw)*std::complex<real_t>(0,1)*(-F[site][0].at(i,j)-F[site][5].at(i,j));

#ifdef ADJOINT
				cloverMatrixInverse[site].at(3*0 + i, 3*1 +j) += (kappa*csw)*std::complex<real_t>(+F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*1 + i, 3*0 +j) += (kappa*csw)*std::complex<real_t>(-F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)-F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*2 + i, 3*3 +j) -= (kappa*csw)*std::complex<real_t>(-F[site][1].at(i,j)+F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*3 + i, 3*2 +j) -= (kappa*csw)*std::complex<real_t>(+F[site][1].at(i,j)-F[site][4].at(i,j),(+F[site][2].at(i,j)+F[site][3].at(i,j)));
#endif
#ifndef ADJOINT
				cloverMatrixInverse[site].at(3*0 + i, 3*1 +j) += (kappa*csw)*(+F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*1 + i, 3*0 +j) += (kappa*csw)*(-F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)-F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*2 + i, 3*3 +j) -= (kappa*csw)*(-F[site][1].at(i,j)+F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)));
				cloverMatrixInverse[site].at(3*3 + i, 3*2 +j) -= (kappa*csw)*(+F[site][1].at(i,j)-F[site][4].at(i,j)+multiply_by_I(+F[site][2].at(i,j)+F[site][3].at(i,j)));
#endif
			}
		}
		
		cloverMatrixInverse[site] = inverse(cloverMatrixInverse[site]);
	}
}

} /* namespace Update */
