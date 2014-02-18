/*
 * StoutSmearing.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#include "StoutSmearing.h"
#include "WilsonGaugeAction.h"

namespace Update {

StoutSmearing::StoutSmearing() { }

StoutSmearing::~StoutSmearing() { }

void StoutSmearing::smearing(const extended_gauge_lattice_t& input, extended_gauge_lattice_t& output, real_t rho) {
	WilsonGaugeAction wga(0.);
	
	GaugeGroup identity;
	for (int i = 0; i < numberColors; ++i) {
		for (int j = 0; j < numberColors; ++j) {
			if (i == j) identity.at(i,j) = 1.;
			else identity.at(i,j) = 0.;
		}
	}

#pragma omp parallel for
	for (int site = 0; site < input.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			GaugeGroup staple = htrans(wga.staple(input,site,mu));
			GaugeGroup omega = rho*staple*htrans(input[site][mu]);
			GaugeGroup iStout = 0.5*(htrans(omega) - omega);

			GaugeGroup toExp = iStout - (trace(iStout)/static_cast<real_t>(numberColors))*identity;
			output[site][mu] = this->exp(-toExp)*input[site][mu];
		}
	}
	output.updateHalo();
}

void StoutSmearing::spatialSmearing(const extended_gauge_lattice_t& initialInput, extended_gauge_lattice_t& output, unsigned int numberLevels, real_t rho) {
	extended_gauge_lattice_t input = initialInput;
	typedef extended_gauge_lattice_t LT;
	for (unsigned int level = 0; level < numberLevels; ++level) {
#pragma omp parallel for
		for (int site = 0; site < input.localsize; ++site) {
			for (unsigned int mu = 0; mu < 3; ++mu) {
				GaugeGroup staple;
				set_to_zero(staple);
				for (int nu = 0; nu < 3; ++nu) {
					if (nu != mu) {
						staple += input[LT::sup(site,mu)][nu]*htrans(input[LT::sup(site,nu)][mu])*htrans(input[site][nu]);
						staple += htrans(input[LT::sup(LT::sdn(site,nu),mu)][nu])*htrans(input[LT::sdn(site,nu)][mu])*input[LT::sdn(site,nu)][nu];
					}
				}
				GaugeGroup omega = rho*staple*htrans(input[site][mu]);
				GaugeGroup iStout = 0.5*(htrans(omega) - omega);
				GaugeGroup identity;
				for (int i = 0; i < numberColors; ++i) {
					for (int j = 0; j < numberColors; ++j) {
						if (i == j) identity.at(i,j) = 1.;
						else identity.at(i,j) = 0.;
					}
				}
				GaugeGroup toExp = iStout - (trace(iStout)/static_cast<real_t>(numberColors))*identity;
				output[site][mu] = this->exp(-toExp)*input[site][mu];
			}
			output[site][3] = input[site][3];
		}
		output.updateHalo();
		input = output;
	}
}

GaugeGroup StoutSmearing::exp(const GaugeGroup& toExp) {
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

	GaugeGroup id;
	for (unsigned int i = 0; i < numberColors; ++i) {
		for (unsigned int j = 0; j < numberColors; ++j) {
			if (i == j) id.at(i,j) = 1.;
			else id.at(i,j) = 0.;
		}
	}

	GaugeGroup result = -f2*toExp*toExp + f1*std::complex<real_t>(0.0,-1.0)*toExp + f0*id;

	return result;
#endif
#if NUMCOLORS != 3
#ifdef EIGEN
	Eigen::ComplexEigenSolver<GaugeGroup> es(toExp);
	GaugeGroup result;
	result.zeros();
	for (int i = 0; i < numberColors; ++i) {
		result.at(i,i) = std::exp(es.eigenvalues()[i]);
	}
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

} /* namespace Update */
