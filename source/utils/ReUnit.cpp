#include "ReUnit.h"

namespace Update {

ReUnit::ReUnit() { }

ReUnit::~ReUnit() { }

void ReUnit::execute(environment_t& environment) {
	double errorvalue = 0.;
	
#pragma omp parallel for reduction(+:errorvalue)
	for (int site = 0; site < environment.gaugeLinkConfiguration.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
#if NUMCOLORS > 2
			errorvalue += abs(1 - abs(det(environment.gaugeLinkConfiguration[site][mu])));
			for (unsigned int i = 0; i < numberColors; ++i) {
				for (unsigned j = 0; j < i; ++j) {
					std::complex<real_t> proj(0.,0.);
					for (unsigned int k = 0; k < numberColors; ++k) proj += conj(environment.gaugeLinkConfiguration[site][mu].at(j,k))*environment.gaugeLinkConfiguration[site][mu].at(i,k);
					for (unsigned int k = 0; k < numberColors; ++k) environment.gaugeLinkConfiguration[site][mu].at(i,k) -= proj*environment.gaugeLinkConfiguration[site][mu].at(j,k);
				}
				std::complex<real_t> normvec(0.,0.);
				for (unsigned int k = 0; k < numberColors; ++k) normvec += conj(environment.gaugeLinkConfiguration[site][mu].at(i,k))*environment.gaugeLinkConfiguration[site][mu].at(i,k);
				for (unsigned int k = 0; k < numberColors; ++k) environment.gaugeLinkConfiguration[site][mu].at(i,k) /= abs(normvec);
			}
			std::complex<real_t> deter = det(environment.gaugeLinkConfiguration[site][mu]);
			for (unsigned int k = 0; k < numberColors; ++k) environment.gaugeLinkConfiguration[site][mu].at(numberColors-1,k) *= std::complex<real_t>(cos(arg(deter)),-sin(arg(deter)));
#endif
#if NUMCOLORS == 2
			real_t a = imag(environment.gaugeLinkConfiguration[site][mu].at(0,0));
			real_t b = imag(environment.gaugeLinkConfiguration[site][mu].at(0,1));
			real_t c = real(environment.gaugeLinkConfiguration[site][mu].at(0,1));
			real_t d = real(environment.gaugeLinkConfiguration[site][mu].at(0,0));
			real_t norm = sqrt(a*a+b*b+c*c+d*d);
			a = a/norm;
			b = b/norm;
			c = c/norm;
			d = d/norm;
			environment.gaugeLinkConfiguration[site][mu].at(0,0) = std::complex<real_t>(d,a);
			environment.gaugeLinkConfiguration[site][mu].at(0,1) = std::complex<real_t>(c,b);
			environment.gaugeLinkConfiguration[site][mu].at(1,0) = std::complex<real_t>(-c,b);
			environment.gaugeLinkConfiguration[site][mu].at(1,1) = std::complex<real_t>(d,-a);
			errorvalue += abs(1-norm);
#endif
		}
	}
	reduceAllSum(errorvalue);
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
	if (isOutputProcess()) std::cout << "ReUnit::Error in unity condition: " << errorvalue << std::endl;
}

double ReUnit::testUnitarity(const extended_gauge_lattice_t& lattice) {
	double errorvalue = 0.;
	
#pragma omp parallel for reduction(+:errorvalue)
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			errorvalue += abs(1 - abs(det(lattice[site][mu])));
		}
	}
	reduceAllSum(errorvalue);
	return errorvalue;
}

} /* namespace Update */
