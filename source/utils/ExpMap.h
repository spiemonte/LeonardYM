/*
 * ExpMap.h
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#ifndef EXPMAP_H
#define EXPMAP_H
#include "MatrixTypedef.h"
#include "ToString.h"
#include "LieGenerators.h"
#include "ConvertLattice.h"

namespace Update {

class ExponentialMap {
public:
	ExponentialMap() {
		for (int i = 0; i < numberColors; ++i) {
			for (int j = 0; j < numberColors; ++j) {
				if (i == j) identity.at(i,j) = 1.;
				else identity.at(i,j) = 0.;
			}
		}
		for (int i = 0; i < numberColors*numberColors-1; ++i) {
			for (int j = 0; j < numberColors*numberColors-1; ++j) {
				if (i == j) adjoint_identity.at(i,j) = 1.;
				else adjoint_identity.at(i,j) = 0.;
			}
		}
	}

	GaugeGroup log(const GaugeGroup& toLog) const {
#ifdef EIGEN
		Eigen::ComplexEigenSolver<GaugeGroup> es(toLog);
		GaugeGroup result;
		result.zeros();
		for (int i = 0; i < numberColors; ++i) {
			result.at(i,i) = std::log(es.eigenvalues()[i]);
		}
		GaugeGroup m = es.eigenvectors();
		return m * result * htrans(m);
#endif
#ifdef ARMADILLO
		FundamentalVector eigval;
		GaugeGroup eigvec;
		arma::eig_gen(eigval, eigvec, toLog);
		GaugeGroup result;
		set_to_zero(result);
		for (int i = 0; i < numberColors; ++i) {
			result.at(i,i) = std::log(eigval[i]);
		}
		return eigvec * result * htrans(eigvec);
#endif
	}

	ForceVector parameters(const GaugeGroup& toLog) const {
		GaugeGroup lg = this->log(toLog);
		ForceVector result;
		for (int c = 0; c < numberColors*numberColors - 1; ++c) result[c] = 2.*(imag(trace(lieGenerator.get(c)*lg)));
		return result;
	}

	void exp(const ForceVector& toExp, GaugeGroup& result) const {
		GaugeGroup M;
		set_to_zero(M);
		for (int c = 0; c < numberColors*numberColors - 1; ++c) M += lieGenerator.get(c)*std::complex<real_t>(0,toExp[c]);
		result = this->exp(M);
	}
	
	
	void exp(ForceVector& toExp, AdjointGroup& result) const {
		GaugeGroup fundamentalExp;
		this->exp(toExp,fundamentalExp);
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::toAdjoint(fundamentalExp,result);
	}
	
	AdjointGroup exp(const AdjointGroup& toExp) const {
#ifdef EIGEN
		typedef Eigen::Matrix< complex, numberColors*numberColors - 1, numberColors*numberColors - 1 > CFermionGroup;
		Eigen::ComplexEigenSolver<CFermionGroup> es( toExp.cast<std::complex<real_t> >() );
		CFermionGroup result;
		result.zeros();
		real_t eigvs_sum = 0.;
		for (int i = 0; i < numberColors*numberColors - 1; ++i) {
			result.at(i,i) = std::exp(es.eigenvalues()[i]);
			eigvs_sum += fabs(real(conj(es.eigenvalues()[i])*es.eigenvalues()[i]));
		}
		if (eigvs_sum < 0.0000000000000001) return adjoint_identity+toExp+(0.5)*toExp*toExp;
		CFermionGroup m = es.eigenvectors();
		m = (m * result * htrans(m));
		AdjointGroup a_result;
		for (int i = 0; i < numberColors*numberColors - 1; ++i) {
			for (int j = 0; j < numberColors*numberColors - 1; ++j) {
				a_result.at(i,j) = real(m.at(i,j));
			}
		}
		return a_result;
#endif
#ifdef ARMADILLO
		std::cout << "ExpMap::Adjoint exponential not defined with armadillo!" << std::endl;
		exit(5);
#endif
	}

	GaugeGroup exp(const GaugeGroup& toExp) const {
#if NUMCOLORS == 3
		real_t c0,c1,c0max,theta,u,w,ksi; //TODO: temporarily double (I've been unable to use Real)
		std::complex<real_t> f0,f1,f2,h0,h1,h2;
		bool inv_c0=0;
		c0=real(det(std::complex<real_t>(0.0,-1.0)*toExp));
		if (c0<0){
			c0=-c0;
			inv_c0=1;
		}
		c1=-0.5*real(trace(toExp*toExp));
		c0max=(2.*c1/3.)*sqrt(c1/3.);
		theta=acos(c0/c0max);
		u=sqrt(c1/3.)*cos(theta/3.);
		w=sqrt(c1)*sin(theta/3.);

		std::complex<real_t> dumexp2(cos(2.*u),sin(2.*u));
		std::complex<real_t> dumexp1(cos(u),-sin(u));

		if(abs(w)<0.05) ksi=1.-(1./6.)*w*w*(1.-(1./20.)*w*w*(1.-(1./42.)*w*w));
		else ksi=sin(w)/w;

		h0=(u*u-w*w)*dumexp2+dumexp1*(8.*u*u*cos(w)+2.*u*std::complex<real_t>(0.0,1.0)*(3.*u*u+w*w)*ksi);
		h1=2.*u*dumexp2-dumexp1*(2.*u*cos(w)+std::complex<real_t>(0.0,-1.0)*(3.*u*u-w*w)*ksi);
		h2=dumexp2-dumexp1*(cos(w)+3.*u*std::complex<real_t>(0.0,1.0)*ksi);

		f0=h0/(9.*u*u-w*w);
		f1=h1/(9.*u*u-w*w);
		f2=h2/(9.*u*u-w*w);
	
		if(inv_c0){
				f0=conj(f0);
				f1=-conj(f1);
				f2=conj(f2);
		}
	

		GaugeGroup result = -f2*toExp*toExp + f1*std::complex<real_t>(0.0,-1.0)*toExp + f0*identity;
	
		return result;
#endif
#if NUMCOLORS == 2
		real_t a = imag(toExp.at(0,0));
		real_t b = real(toExp.at(0,1));
		real_t c = imag(toExp.at(0,1));
		real_t rho = sqrt(a*a+b*b+c*c);
		std::complex<real_t> phi(cos(rho),sin(rho));
	
		GaugeGroup res;
		res.at(0,0) = (-(a-rho)*conj(phi)+(a+rho)*phi)/(2.*rho);
		res.at(0,1) = (std::complex<real_t>(-c,b)*conj(phi)+std::complex<real_t>(c,-b)*phi)/(2.*rho);
		res.at(1,0) = -conj(res.at(0,1));
		res.at(1,1) = conj(res.at(0,0));
		return res;
#endif
#if NUMCOLORS != 3 && NUMCOLORS != 2
#ifdef EIGEN
		Eigen::ComplexEigenSolver<GaugeGroup> es(toExp);
		GaugeGroup result;
		result.zeros();
		real_t eigvs_sum = 0.;
		for (int i = 0; i < numberColors; ++i) {
			result.at(i,i) = std::exp(es.eigenvalues()[i]);
			eigvs_sum += fabs(real(conj(es.eigenvalues()[i])*es.eigenvalues()[i]));
		}
		if (eigvs_sum < 0.0000000000000001) return identity+toExp+(0.5)*toExp*toExp;
		GaugeGroup m = es.eigenvectors();
		return m * result * htrans(m);
#endif
#ifdef ARMADILLO
		FundamentalVector eigval;
		GaugeGroup eigvec;
		arma::eig_gen(eigval, eigvec, toExp);
		GaugeGroup result;
		set_to_zero(result);
		for (int i = 0; i < numberColors; ++i) {
			result.at(i,i) = std::exp(eigval[i]);
		}
		return eigvec * result * htrans(eigvec);
#endif
#endif
	}

protected:
	GaugeGroup identity;
	AdjointGroup adjoint_identity;
	LieGenerator<GaugeGroup> lieGenerator;
	LieGenerator<AdjointGroup> adjointLieGenerator;
};

} /* namespace Update */

#endif
