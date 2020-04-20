#include "BasicDiracWilsonOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"

namespace Update {

BasicDiracWilsonOperator::BasicDiracWilsonOperator() : DiracOperator() { }

BasicDiracWilsonOperator::BasicDiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

BasicDiracWilsonOperator::~BasicDiracWilsonOperator() { }

void BasicDiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {		
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

		
		if (gamma5) {
			//The final result is gamma5*input - kappa*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = std::complex<real_t>(real(input[site][0][n])-kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = std::complex<real_t>(real(input[site][1][n])-kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = std::complex<real_t>(-(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n]))),-(imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n]))));
				output[site][3][n] = std::complex<real_t>(-(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n]))),-(imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n]))));
			}
		}
		else {
			//The final result is input - kappa*gamma5*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = std::complex<real_t>(real(input[site][0][n]) - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(input[site][0][n]) - kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = std::complex<real_t>(real(input[site][1][n]) - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(input[site][1][n]) - kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = std::complex<real_t>(real(input[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])),imag(input[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
				output[site][3][n] = std::complex<real_t>(real(input[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])),imag(input[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
			}
		}
	}
	output.updateHalo();
}

void BasicDiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	typedef reduced_fermion_lattice_t Lattice;
	typedef reduced_dirac_vector_t Vector;

#pragma omp parallel for
	for (int site = 0; site < lattice.localsize; ++site) {
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

		
		if (gamma5) {
			//The final result is gamma5*vector1 - kappa*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(real(vector1[site][0][n])-kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(vector1[site][0][n])-kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(real(vector1[site][1][n])-kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(vector1[site][1][n])-kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = alpha*vector2[site][2][n] + std::complex<real_t>(-(real(vector1[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n]))),-(imag(vector1[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n]))));
				output[site][3][n] = alpha*vector2[site][3][n] + std::complex<real_t>(-(real(vector1[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n]))),-(imag(vector1[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n]))));
			}
		}
		else {
			//The final result is vector1 - kappa*gamma5*hopping
			for (int n = 0; n < diracVectorLength; ++n) {
				output[site][0][n] = alpha*vector2[site][0][n] + std::complex<real_t>(real(vector1[site][0][n]) - kappa*(real(tmp_plus[0][0][n])+real(tmp_plus[1][0][n])+real(tmp_plus[2][0][n])+real(tmp_plus[3][0][n])),imag(vector1[site][0][n]) - kappa*(imag(tmp_plus[0][0][n])+imag(tmp_plus[1][0][n])+imag(tmp_plus[2][0][n])+imag(tmp_plus[3][0][n])));
				output[site][1][n] = alpha*vector2[site][1][n] + std::complex<real_t>(real(vector1[site][1][n]) - kappa*(real(tmp_plus[0][1][n])+real(tmp_plus[1][1][n])+real(tmp_plus[2][1][n])+real(tmp_plus[3][1][n])),imag(vector1[site][1][n]) - kappa*(imag(tmp_plus[0][1][n])+imag(tmp_plus[1][1][n])+imag(tmp_plus[2][1][n])+imag(tmp_plus[3][1][n])));
				output[site][2][n] = alpha*vector2[site][2][n] + std::complex<real_t>(real(vector1[site][2][n]) + kappa*(real(tmp_minus[1][1][n]) + real(tmp_minus[3][0][n]) - imag(tmp_minus[0][1][n]) - imag(tmp_minus[2][0][n])),imag(vector1[site][2][n]) + kappa*(real(tmp_minus[0][1][n]) + real(tmp_minus[2][0][n]) + imag(tmp_minus[1][1][n]) + imag(tmp_minus[3][0][n])));
				output[site][3][n] = alpha*vector2[site][3][n] + std::complex<real_t>(real(vector1[site][3][n]) + kappa*(imag(tmp_minus[2][1][n]) - imag(tmp_minus[0][0][n]) - real(tmp_minus[1][0][n]) + real(tmp_minus[3][1][n])),imag(vector1[site][3][n]) + kappa*(real(tmp_minus[0][0][n]) - real(tmp_minus[2][1][n]) - imag(tmp_minus[1][0][n]) + imag(tmp_minus[3][1][n])));
			}
		}
	}
	output.updateHalo();//TODO is needed?
}

FermionForce* BasicDiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

} /* namespace Update */
