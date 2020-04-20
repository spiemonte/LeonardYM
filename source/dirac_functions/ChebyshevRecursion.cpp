#include "ChebyshevRecursion.h"

namespace Update {

ChebyshevRecursion::ChebyshevRecursion() : a(0), b(1), n(2) { }

ChebyshevRecursion::ChebyshevRecursion(real_t _a, real_t _b, unsigned int _n) : a(_a), b(_b), n(_n) { }

ChebyshevRecursion::ChebyshevRecursion(const ChebyshevRecursion& copy) : a(copy.a), b(copy.b), n(copy.n) { }

ChebyshevRecursion::~ChebyshevRecursion() { }

void ChebyshevRecursion::evaluate(DiracOperator* diracOperator, reduced_dirac_vector_t& output, const reduced_dirac_vector_t& input) {
	reduced_dirac_vector_t *next = new reduced_dirac_vector_t();
	reduced_dirac_vector_t *previous =  new reduced_dirac_vector_t(input);
	reduced_dirac_vector_t *actual =  new reduced_dirac_vector_t();
	reduced_dirac_vector_t *swap;

	diracOperator->multiply(*actual, input);
#pragma omp parallel for
	for (int site = 0; site < actual->completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			(*actual)[site][mu] = (2./(b-a))*(*actual)[site][mu] - ((a+b)/(b-a))*input[site][mu];
		}
	}

	for (unsigned int m = 0; m < n; ++m) {
		//*next = 2.*(2.*input/(b-a) - (a+b)/(b-a))*(*actual) - *previous;
		diracOperator->multiply(*next, *actual);

#pragma omp parallel for
		for (int site = 0; site < next->completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*next)[site][mu] = 2.*((2./(b-a))*(*next)[site][mu] - ((a+b)/(b-a))*(*actual)[site][mu]) - (*previous)[site][mu];
			}
		}

		swap = previous;
		previous = actual;
		actual = next;
		next = swap;

	}
	
	output = (*actual);

	delete previous;
	delete next;
	delete actual;
/*
	real_t sigma = a/(nu-b);
	diracOperator->multiplyAdd(*actual, input, input, -b);
#pragma omp parallel for
	for (int site = 0; site < next->completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			(*actual)[site][mu] = (sigma/a)*(*actual)[site][mu];
		}
	}
	
	for (unsigned int m = 0; m < n; ++m) {
		real_t sigma_next = 1./(2.*(nu-b)/a - sigma);
		diracOperator->multiplyAdd(*next, *actual, *actual, -b);
		
#pragma omp parallel for
		for (int site = 0; site < next->completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				(*next)[site][mu] = (2.*sigma_next/a)*(*next)[site][mu] - sigma*sigma_next*(*previous)[site][mu];
			}
		}
		

		swap = previous;
		previous = actual;
		actual = next;
		next = swap;
		sigma = sigma_next;
	}

	output = (*actual);

	delete previous;
	delete next;
	delete actual;*/
}



complex ChebyshevRecursion::evaluate(const complex& input) const {
	complex *next = new complex();
	complex *previous = new complex(input);
	complex *actual = new complex();
	complex *swap;

	*actual = 2.*input/(b-a) - (a+b)/(b-a);
	*previous = 1.;
	for (unsigned int m = 0; m < n; ++m) {
		*next = 2.*(2.*input/(b-a) - (a+b)/(b-a))*(*actual) - *previous;

		swap = previous;
		previous = actual;
		actual = next;
		next = swap;
	}

	complex output = (*actual);

	delete previous;
	delete next;
	delete actual;
	return output;
}

void ChebyshevRecursion::setParameters(real_t _a, real_t _b, unsigned int _n) {
	a = _a;
	b = _b;
	n = _n;
}



} /* namespace Update */
