/*
 * LieGenerators.h
 *
 *  Created on: Nov 28, 2012
 *      Author: spiem_01
 */

#ifndef LIEGENERATORS_H_
#define LIEGENERATORS_H_
#include "MatrixTypedef.h"
#define I std::complex<real_t>(0,1)

namespace Update {

template<typename Representation> class LieGenerator {
public:
	unsigned int numberGenerators() const { return 0; } //TODO
	const GaugeGroup& get(unsigned int c) const { return GaugeGroup(); }//TODO
};

#if NUMCOLORS == 2

template<> class LieGenerator<GaugeGroup> {
public:
	LieGenerator() {
#ifdef EIGEN
		gen[0] << 0, 1, 1, 0;
		gen[1] << 0, -I, I, 0;
		gen[2] << 1, 0, 0, -1;
		gen[0] = gen[0]/2.;
		gen[1] = gen[1]/2.;
		gen[2] = gen[2]/2.;
#endif
#ifdef ARMADILLO
		gen[0] << 0 <<  1 << arma::endr
			   << 1 <<  0;
		gen[1] << 0 << -I << arma::endr
			   << I <<  0;
		gen[2] << 1 <<  0 << arma::endr
			   << 0 << -1;
		gen[0] = gen[0]/2.;
		gen[1] = gen[1]/2.;
		gen[2] = gen[2]/2.;
#endif
#ifdef MATRIXTOOLKIT
		gen[0][0][0] = 0.;
		gen[0][0][1] = 1.;
		gen[0][1][0] = 1.;
		gen[0][1][1] = 0.;

		gen[1][0][0] = 0.;
		gen[1][0][1] = -I;
		gen[1][1][0] = I;
		gen[1][1][1] = 0.;

		gen[2][0][0] = 1.;
		gen[2][0][1] = 0.;
		gen[2][1][0] = 0.;
		gen[2][1][1] = -1.;

		gen[0] = gen[0]/2.;
		gen[1] = gen[1]/2.;
		gen[2] = gen[2]/2.;
#endif
	}

	const GaugeGroup& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 3;
	}
private:
	GaugeGroup gen[3];
};

template<> class LieGenerator<AdjointGroup> {
public:
	LieGenerator() {
#ifdef EIGEN
		gen[0] << 0, 0, 0, 0, 0, -I, 0, I, 0;
		gen[1] << 0, 0, I, 0, 0, 0, -I, 0, 0;
		gen[2] << 0, -I, 0, I, 0, 0, 0, 0, 0;
#endif
#ifdef ARMADILLO
		gen[0] << 0 << 0 <<  0 << arma::endr
			   << 0 << 0 << -I << arma::endr
			   << 0 << I << 0;
		gen[1] <<  0 << 0 << I << arma::endr
			   <<  0 << 0 << 0 << arma::endr
			   << -I << 0 << 0;
		gen[2] << 0 << -I << 0 << arma::endr
			   << I <<  0 << 0 << arma::endr
			   << 0 <<  0 << 0;
#endif
#ifdef MATRIXTOOLKIT
		set_to_zero(gen[0]);
		set_to_zero(gen[1]);
		set_to_zero(gen[2]);

		gen[0][1][2] = -I;
		gen[0][2][1] = I;

		gen[1][0][2] = I;
		gen[1][2][0] = -I;

		gen[2][0][1] = -I;
		gen[2][1][0] = I;
#endif
	}

	const FermionicForceMatrix& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 3;
	}
private:
	FermionicForceMatrix gen[3];
};

#endif

#if NUMCOLORS == 3

template<> class LieGenerator<GaugeGroup> {
public:
	LieGenerator() {
		gen[0] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
		gen[1] << 0, 0, 1, 0, 0, 0, 1, 0, 0;
		gen[2] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
		gen[3] << 0, -I, 0, I, 0, 0, 0, 0, 0;
		gen[4] << 0, 0, -I, 0, 0, 0, I, 0, 0;
		gen[5] << 0, 0, 0, 0, 0, -I, 0, I, 0;
		gen[6] << 1, 0, 0, 0, -1, 0, 0, 0, 0;
		gen[7] << 1./(sqrt(3.)), 0, 0, 0, 1./(sqrt(3.)), 0, 0, 0, -2./(sqrt(3.));
		for (unsigned int i = 0; i < 8; ++i) gen[i] = gen[i]/2.;
	}

	const GaugeGroup& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 8;
	}
private:
	GaugeGroup gen[8];
};

template<> class LieGenerator<AdjointGroup> {
public:
	LieGenerator() {
		gen[0] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, -I, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, I, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[1] << 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), -((I*sqrt(3.))/2.), -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, (I*sqrt(3.))/2., 0, 0, 0;
		gen[2] << 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, I/2., -((I*sqrt(3.))/2.), 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, (I*sqrt(3.))/2., 0, 0;
		gen[3] << 0, 0, 0, 0, 0, 0, I, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, -I, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[4] << 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, I/2., (I*sqrt(3.))/2., I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, -((I*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0;
		gen[5] << 0, -(I/2.), 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), (I*sqrt(3.))/2., 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, -((I*sqrt(3.))/2.), 0, 0, 0, 0, 0;
		gen[6] << 0, 0, 0, -I, 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, I/2., 0, 0, I, 0, 0, 0, 0, 0, 0, 0, 0, I/2., 0, 0, 0, 0, 0, 0, 0, 0, -(I/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[7] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -((I*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0, 0, 0, -((I*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (I*sqrt(3.))/2., 0, 0, 0, 0, 0, 0, 0, 0, (I*sqrt(3.))/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	}

	const FermionicForceMatrix& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 8;
	}
private:
	FermionicForceMatrix gen[8];
};

#endif


}


#endif /* LIEGENERATORS_H_ */
