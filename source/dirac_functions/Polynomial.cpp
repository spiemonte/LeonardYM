/*
 * Polynomial.cpp
 *
 *  Created on: Apr 30, 2012
 *      Author: spiem_01
 */

#include "Polynomial.h"

namespace Update {

Polynomial::Polynomial() { }

Polynomial::Polynomial(const std::vector< complex >& _roots, const complex& _scaling) : roots(_roots), scaling(_scaling) { }

Polynomial::~Polynomial() { }

void Polynomial::evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& original_output, const extended_dirac_vector_t& original_input) {
	//We work with reduced halos
	reduced_dirac_vector_t output;
	reduced_dirac_vector_t input = original_input;
	
	if (roots.size() % 2 == 0) {//TODO minus sign
		std::vector<complex>::iterator i = roots.begin();

		//First step
		diracOperator->multiplyAdd(tmp1, input, input, -*i);
#pragma omp parallel for
		for (int site = 0; site < tmp1.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = fabs(real(scaling))*tmp1[site][mu];
		}
		
		++i;
		
		while (i != roots.end()) {
			diracOperator->multiplyAdd(tmp2, tmp1, tmp1, -*i);
#pragma omp parallel for
			for (int site = 0; site < tmp2.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) tmp2[site][mu] = fabs(real(scaling))*tmp2[site][mu];
			}
			
			++i;
			if (i != roots.end()) {
				diracOperator->multiplyAdd(tmp1, tmp2, tmp2, -*i);
#pragma omp parallel for
				for (int site = 0; site < tmp1.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = fabs(real(scaling))*tmp1[site][mu];
				}
				
				++i;
			}
		}
		
		if (real(scaling) < 0) {
			for (int site = 0; site < input.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] = -tmp2[site][mu];
				}
			}
			output.updateHalo();
		}
		else {
			output = tmp2;
		}
	} else {
		std::vector<complex>::iterator i = roots.begin();

		diracOperator->multiplyAdd(tmp1, input, input, -*i);
#pragma omp parallel for
		for (int site = 0; site < tmp1.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = tmp1[site][mu]*fabs(real(scaling));
		}
		
		++i;
		while (i != roots.end()) {
			diracOperator->multiplyAdd(tmp2, tmp1, tmp1, -*i);
#pragma omp parallel for
			for (int site = 0; site < tmp2.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) tmp2[site][mu] = tmp2[site][mu]*fabs(real(scaling));
			}
			
			++i;
			if (i != roots.end()) {
				diracOperator->multiplyAdd(tmp1, tmp2, tmp2, -*i);
#pragma omp parallel for
				for (int site = 0; site < tmp1.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = tmp1[site][mu]*fabs(real(scaling));
				}
				
				++i;
			}
		}

		if (real(scaling) < 0) {
			for (int site = 0; site < input.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] = -tmp1[site][mu];
				}
			}
			output.updateHalo();
		}
		else {
			output = tmp1;
		}
	}
	original_output = output;
}

void Polynomial::evaluate(DiracOperator* diracOperator, extended_dirac_vector_t& original_output, const extended_dirac_vector_t& original_input, DiracOperator* preconditioner) {
	//We work with reduced halos
	reduced_dirac_vector_t output;
	reduced_dirac_vector_t input = original_input;
	
	if (roots.size() % 2 == 0) {//TODO minus sign
		std::vector<complex>::iterator i = roots.begin();

		//First step
		diracOperator->multiply(tmp2, input);
		preconditioner->multiplyAdd(tmp1, tmp2, input, -*i);
#pragma omp parallel for
		for (int site = 0; site < tmp1.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = fabs(real(scaling))*tmp1[site][mu];
		}
		
		++i;
		
		while (i != roots.end()) {
			diracOperator->multiply(output, tmp1);
			preconditioner->multiplyAdd(tmp2, output, tmp1, -*i);
#pragma omp parallel for
			for (int site = 0; site < tmp2.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) tmp2[site][mu] = fabs(real(scaling))*tmp2[site][mu];
			}
			
			++i;
			if (i != roots.end()) {
				diracOperator->multiply(output, tmp2);
				preconditioner->multiplyAdd(tmp1, output, tmp2, -*i);
#pragma omp parallel for
				for (int site = 0; site < tmp1.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = fabs(real(scaling))*tmp1[site][mu];
				}
				
				++i;
			}
		}
		
		if (real(scaling) < 0) {
			preconditioner->multiply(output, tmp2);
			for (int site = 0; site < input.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] = -output[site][mu];
				}
			}
		}
		else {
			preconditioner->multiply(output, tmp2);
		}
	} else {
		std::vector<complex>::iterator i = roots.begin();

		diracOperator->multiply(tmp2, input);
		preconditioner->multiplyAdd(tmp1, tmp2, input, -*i);
#pragma omp parallel for
		for (int site = 0; site < tmp1.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = tmp1[site][mu]*fabs(real(scaling));
		}
		
		++i;
		while (i != roots.end()) {
			diracOperator->multiply(output, tmp1);
			preconditioner->multiplyAdd(tmp2, output, tmp1, -*i);
#pragma omp parallel for
			for (int site = 0; site < tmp2.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) tmp2[site][mu] = tmp2[site][mu]*fabs(real(scaling));
			}
			
			++i;
			if (i != roots.end()) {
				diracOperator->multiply(output, tmp2);
				preconditioner->multiplyAdd(tmp1, output, tmp2, -*i);
#pragma omp parallel for
				for (int site = 0; site < tmp1.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) tmp1[site][mu] = tmp1[site][mu]*fabs(real(scaling));
				}
				
				++i;
			}
		}

		if (real(scaling) < 0) {
			preconditioner->multiply(output, tmp1);
			for (int site = 0; site < input.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					output[site][mu] = -output[site][mu];
				}
			}
		}
		else {
			preconditioner->multiply(output, tmp1);
		}
	}
	original_output = output;
}

complex Polynomial::evaluate(const complex& x) const {
	complex result;
	if (roots.size() % 2 == 0) {
		std::vector<complex>::const_iterator i = roots.begin();
		result = static_cast<real_t>(fabs(real(scaling)))*(x - *i);
		++i;
		while (i != roots.end()) {
			result = static_cast<real_t>(fabs(real(scaling)))*result*(x - *i);
			++i;
			if (i != roots.end()) {
				result = static_cast<real_t>(fabs(real(scaling)))*result*(x - *i);
				++i;
			}
		}
		if (real(scaling) < 0) return -result;
		else return result;
	} else {
		std::vector<complex>::const_iterator i = roots.begin();
		result = static_cast<real_t>(fabs(real(scaling)))*(x - *i);
		++i;
		while (i != roots.end()) {
			result = static_cast<real_t>(fabs(real(scaling)))*result*(x - *i);
			++i;
			if (i != roots.end()) {
				result = static_cast<real_t>(fabs(real(scaling)))*result*(x - *i);
				++i;
			}
		}
		if (real(scaling) < 0) return -result;
		else return result;
	}
}

void Polynomial::setScaling(const complex& _scaling) {
	scaling = _scaling;
}

complex Polynomial::getScaling() const {
	return scaling;
}

void Polynomial::setRoots(const std::vector<complex>& _roots) {
	roots = _roots;
}

std::vector<complex> Polynomial::getRoots() const {
	return roots;
}

} /* namespace Update */
