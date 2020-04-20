#include "DiracWilsonOperator.h"
#include "hmc_forces/DiracWilsonFermionForce.h"

namespace Update {

DiracWilsonOperator::DiracWilsonOperator() : DiracOperator() { }

DiracWilsonOperator::DiracWilsonOperator(const extended_fermion_lattice_t& _lattice, real_t _kappa, bool _gamma5) : DiracOperator(_lattice, _kappa, _gamma5) { }

DiracWilsonOperator::~DiracWilsonOperator() { }

inline real_t conj(const real_t& t) {
	return t;
}

void DiracWilsonOperator::multiply(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
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

 void DiracWilsonOperator::multiplyAdd(reduced_dirac_vector_t& output, const reduced_dirac_vector_t& vector1, const reduced_dirac_vector_t& vector2, const complex& alpha) {
	 typedef reduced_fermion_lattice_t Lattice;
	 typedef reduced_dirac_vector_t Vector;

	 const reduced_fermion_lattice_t& linkconf = (lattice);

#pragma omp parallel for
	 for (int site = 0; site < lattice.localsize; ++site) {
		 //the most efficient unrolling of the dirac operator
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

FermionForce* DiracWilsonOperator::getForce() const {
	return new DiracWilsonFermionForce(kappa);
}

} /* namespace Update */
