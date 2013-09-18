/*
 * Gamma.h
 *
 *  Created on: Nov 30, 2012
 *      Author: spiem_01
 */

#ifndef GAMMA_H_
#define GAMMA_H_
#include "MatrixTypedef.h"

namespace Update {

class Gamma {
public:
	//This function returns the value of (\gamma_5(1+\gamma_\mu))_{\alpha,\beta}
	static std::complex<real_t> g5idpg(unsigned int mu, unsigned int alpha, unsigned int beta) {
		return gamma5_id_plus_gamma[mu][alpha][beta];
	}

	//This function returns the value of (\gamma_5(1-\gamma_\mu))_{\alpha,\beta}
	static std::complex<real_t> g5idmg(unsigned int mu, unsigned int alpha, unsigned int beta) {
		return gamma5_id_minus_gamma[mu][alpha][beta];
	}

	static std::complex<real_t> gamma5(unsigned int alpha, unsigned int beta) {
		return _gamma5[alpha][beta];
	}
private:
	const static std::complex<real_t> gamma5_id_plus_gamma[4][4][4];
	const static std::complex<real_t> gamma5_id_minus_gamma[4][4][4];
	const static std::complex<real_t> _gamma5[4][4];
};

class Sigma {
public:
	//This function returns the value of (\gamma_5(\sigma_{\mu\nu}))_{\alpha,\beta}
	static std::complex<real_t> g5sigma(unsigned int mu, unsigned int nu, unsigned int alpha, unsigned int beta) {
		return gamma5_sigma[mu][nu][alpha][beta];
	}

	static std::complex<real_t> sigma(unsigned int mu, unsigned int nu, unsigned int alpha, unsigned int beta) {
		return _sigma[mu][nu][alpha][beta];
	}

private:
	const static std::complex<real_t> gamma5_sigma[4][4][4][4];
	const static std::complex<real_t> _sigma[4][4][4][4];
};

}


#endif /* GAMMA_H_ */
