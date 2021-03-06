#ifndef LIEGENERATORS_H_
#define LIEGENERATORS_H_
#include "MatrixTypedef.h"
#define II std::complex<real_t>(0,1)

namespace Update {

template<typename Representation> class LieGenerator {
public:
	virtual unsigned int numberGenerators() const = 0;// { return 0; } //TODO
	virtual const GaugeGroup& get(unsigned int c) const = 0;//{ return GaugeGroup(); }//TODO
};

#if NUMCOLORS == 2

template<> class LieGenerator<GaugeGroup> {
public:
	LieGenerator() {
#ifdef EIGEN
		gen[0] << 0, 1, 1, 0;
		gen[1] << 0, -II, II, 0;
		gen[2] << 1, 0, 0, -1;
		gen[0] = gen[0]/2.;
		gen[1] = gen[1]/2.;
		gen[2] = gen[2]/2.;
#endif
#ifdef ARMADILLO
		gen[0] << 0 <<  1 << arma::endr
			   << 1 <<  0;
		gen[1] << 0 << -II << arma::endr
			   << II <<  0;
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
		gen[1][0][1] = -II;
		gen[1][1][0] = II;
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
		gen[0] << 0, 0, 0, 0, 0, -II, 0, II, 0;
		gen[1] << 0, 0, II, 0, 0, 0, -II, 0, 0;
		gen[2] << 0, -II, 0, II, 0, 0, 0, 0, 0;
#endif
#ifdef ARMADILLO
		gen[0] << 0 << 0 <<  0 << arma::endr
			   << 0 << 0 << -II << arma::endr
			   << 0 << II << 0;
		gen[1] <<  0 << 0 << II << arma::endr
			   <<  0 << 0 << 0 << arma::endr
			   << -II << 0 << 0;
		gen[2] << 0 << -II << 0 << arma::endr
			   << II <<  0 << 0 << arma::endr
			   << 0 <<  0 << 0;
#endif
#ifdef MATRIXTOOLKIT
		set_to_zero(gen[0]);
		set_to_zero(gen[1]);
		set_to_zero(gen[2]);

		gen[0][1][2] = -II;
		gen[0][2][1] = II;

		gen[1][0][2] = II;
		gen[1][2][0] = -II;

		gen[2][0][1] = -II;
		gen[2][1][0] = II;
#endif
	}

	const Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 >& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 3;
	}
private:
	Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 > gen[3];
};

//#endif

#elif NUMCOLORS == 3

template<> class LieGenerator<GaugeGroup> {
public:
	LieGenerator() {
		gen[0] << 0, 1, 0, 1, 0, 0, 0, 0, 0;
		gen[1] << 0, 0, 1, 0, 0, 0, 1, 0, 0;
		gen[2] << 0, 0, 0, 0, 0, 1, 0, 1, 0;
		gen[3] << 0, -II, 0, II, 0, 0, 0, 0, 0;
		gen[4] << 0, 0, -II, 0, 0, 0, II, 0, 0;
		gen[5] << 0, 0, 0, 0, 0, -II, 0, II, 0;
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
		gen[0] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, -II, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, II, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[1] << 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), -((II*sqrt(3.))/2.), -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, (II*sqrt(3.))/2., 0, 0, 0;
		gen[2] << 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, II/2., -((II*sqrt(3.))/2.), 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, (II*sqrt(3.))/2., 0, 0;
		gen[3] << 0, 0, 0, 0, 0, 0, II, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, -II, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[4] << 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, II/2., (II*sqrt(3.))/2., II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, -((II*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0;
		gen[5] << 0, -(II/2.), 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), (II*sqrt(3.))/2., 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, -((II*sqrt(3.))/2.), 0, 0, 0, 0, 0;
		gen[6] << 0, 0, 0, -II, 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, II/2., 0, 0, II, 0, 0, 0, 0, 0, 0, 0, 0, II/2., 0, 0, 0, 0, 0, 0, 0, 0, -(II/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
		gen[7] << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -((II*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0, 0, 0, -((II*sqrt(3.))/2.), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (II*sqrt(3.))/2., 0, 0, 0, 0, 0, 0, 0, 0, (II*sqrt(3.))/2., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
	}

	const Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 >& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return 8;
	}
private:
	Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 > gen[8];
};

#else

template<> class LieGenerator<GaugeGroup> {
public:
	LieGenerator() {
		int index = 0;
		for (int i = 0; i < numberColors; ++i) {
			for (int j = 0; j < i; ++j) {
				for (int m = 0; m < numberColors; ++m) {
					for (int n = 0; n < numberColors; ++n) {
						if (i == m && j == n) gen[index].at(m,n) = 1;
						else gen[index].at(m,n) = 0;
						if (i == n && j == m) gen[index].at(m,n) += 1;
					}
				}
				++index;
			}
		}
		for (int i = 0; i < numberColors; ++i) {
			for (int j = 0; j < i; ++j) {
				for (int m = 0; m < numberColors; ++m) {
					for (int n = 0; n < numberColors; ++n) {
						if (i == m && j == n) gen[index].at(m,n) = II;
						else gen[index].at(m,n) = 0;
						if (i == n && j == m) gen[index].at(m,n) -= II;
					}
				}
				++index;
			}
		}
		for (int i = 1; i < numberColors; ++i) {
			for (int m = 0; m < numberColors; ++m) {
				for (int n = 0; n < numberColors; ++n) {
					if (m == n && m <= i) gen[index].at(m,n) = sqrt(2./static_cast<real_t>((i+1)*i));
					else gen[index].at(m,n) = 0;
					if (m == n && m == i) gen[index].at(m,n) -= (i+1)*sqrt(2./static_cast<real_t>((i+1)*i));
				}
			}
			++index;
		}
		for (unsigned int i = 0; i < numberColors*numberColors - 1; ++i) gen[i] = gen[i]/2.;
	}

	const GaugeGroup& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return numberColors*numberColors-1;
	}
private:
	GaugeGroup gen[numberColors*numberColors-1];
};

template<> class LieGenerator<AdjointGroup> {
public:
	LieGenerator() {
		LieGenerator<GaugeGroup> fundamentalGenerators;
		for (unsigned int i = 0; i < numberColors*numberColors - 1; ++i) {
			for (int m = 0; m < numberColors*numberColors - 1; ++m) {
				for (int n = 0; n < numberColors*numberColors - 1; ++n) {
					if (m != n) gen[i].at(m,n) = -2.*trace((fundamentalGenerators.get(m)*fundamentalGenerators.get(n) - fundamentalGenerators.get(n)*fundamentalGenerators.get(m))*fundamentalGenerators.get(i));
					else gen[i].at(m,n) = 0;
				}
			}
		}
	}

	const Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 >& get(unsigned int c) const {
		return gen[c];
	}

	unsigned int numberGenerators() const {
		return numberColors*numberColors-1;
	}
private:
	Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 > gen[numberColors*numberColors-1];
};

#endif

}


#endif /* LIEGENERATORS_H_ */
