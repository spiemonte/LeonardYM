/*
 * MMMRMultishiftSolver.cpp
 *
 *  Created on: Oct 30, 2012
 *      Author: spiem_01
 */

#include "MMMRMultishiftSolver.h"
#include "AlgebraUtils.h"
#include "BiConjugateGradient.h"
#ifdef TEST_PAPI_SPEED
#include <papi.h>


static void test_fail(char *file, int line, char *call, int retval) {
	printf("%s\tFAILED\nLine # %d\n", file, line);

    if ( retval == PAPI_ESYS ) {
        char buf[128];
        memset( buf, '\0', sizeof(buf) );
        sprintf(buf, "System error in %s:", call );
        perror(buf);
    }

    else if ( retval > 0 ) {
        printf("Error calculating: %s\n", call );
    }

    else {
        char errstring[PAPI_MAX_STR_LEN];
        //PAPI_perror(retval, errstring, PAPI_MAX_STR_LEN );
        std::cout << "Papi fatal error! " << std::endl;
        //printf("Error in %s: %s\n", call, errstring );
    }
    printf("\n");
    exit(1);
}

#endif

namespace Update {

MMMRMultishiftSolver::MMMRMultishiftSolver(real_t _epsilon, unsigned int _maxSteps) : MultishiftSolver(_epsilon, _maxSteps), omega(0.95) { }

MMMRMultishiftSolver::~MMMRMultishiftSolver() { }

bool MMMRMultishiftSolver::solve(DiracOperator* dirac, const extended_dirac_vector_t& original_source, std::vector<extended_dirac_vector_t>& original_solutions, const std::vector<real_t>& shifts) {
	//We work with reduced halos
	reduced_dirac_vector_t source = original_source;
	std::vector<reduced_dirac_vector_t> solutions(original_solutions.size());

#ifdef TEST_PAPI_SPEED
	float real_time, proc_time, mflops;
	long long flpins;
	int retval;
	if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
#endif
	
	//The index of the first solution
	unsigned int firstSolution = 0;
	unsigned int length = shifts.size();
	bool flag = true;
	//consider first the negative shifts
	while (shifts[firstSolution] < 0) {
		//Solve it separately
		///
		reduced_dirac_vector_t tmp;
		
#pragma omp parallel for
		for (int site = 0; site < r.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				set_to_zero(solutions[firstSolution][site][mu]);
				p[site][mu] = source[site][mu];
				r[site][mu] = source[site][mu];
			}
		}
		solutions[firstSolution].updateHalo();
		p.updateHalo();
		r.updateHalo();

		long_real_t normResidual = AlgebraUtils::squaredNorm(r);
		for (unsigned int i = 0; i < maxSteps; ++i) {
			dirac->multiplyAdd(tmp, p, p, shifts[firstSolution]);
			long_real_t alpha = normResidual/real(AlgebraUtils::dot(p,tmp));

#pragma omp parallel for
			for (int site = 0; site < r.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solutions[firstSolution][site][mu] = solutions[firstSolution][site][mu] + alpha*p[site][mu];
					r[site][mu] = r[site][mu] - alpha*tmp[site][mu];
				}
			}
			solutions[firstSolution].updateHalo();
			r.updateHalo();//TODO maybe not needed

			long_real_t error = AlgebraUtils::squaredNorm(r);

			if (error < epsilon) {
#ifdef MULTISHIFTLOG
				std::cout << "Multishift BiCGStab steps: " << i << " for shift " << *shift << std::endl;
#endif
				break;
			}

			long_real_t beta = error / normResidual;

#pragma omp parallel for
			for (int site = 0; site < r.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					p[site][mu] = r[site][mu] + beta*p[site][mu];
				}
			}
			p.updateHalo();

			normResidual = error;

			if (i == maxSteps - 1) flag = false;
		}

		///
		++firstSolution;
	}

	std::complex<real_t>* f = new std::complex<real_t>[shifts.size()];
	for (unsigned int i = 0; i < shifts.size(); ++i) f[i] = 1.;

	//Set the solutions to zero
	for (unsigned int i = firstSolution; i < shifts.size(); ++i) {
		AlgebraUtils::setToZero(solutions[i]);
	}

	//Set r to b
#pragma omp parallel for
	for (int site = 0; site < r.completesize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			r[site][mu] = source[site][mu];
		}
	}
	//r.updateHalo(); not needed any more

	for (unsigned int i = 0; i < maxSteps; ++i) {
		//Multiply by the last shifts
		dirac->multiplyAdd(p,r,r,shifts.back());
		std::complex<real_t> alpha = omega*static_cast< std::complex<real_t> >((AlgebraUtils::dot(p,r))/(AlgebraUtils::dot(p,p)));

		//Set the last solution
#pragma omp parallel for
		for (int site = 0; site < r.localsize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solutions[length - 1][site][mu] = solutions[length - 1][site][mu] + alpha*r[site][mu];
			}
		}
		//solutions[length - 1].updateHalo();//TODO maybe not needed

		for (unsigned int index = firstSolution; index < shifts.size() - 1; ++index) {
			f[index] = f[index]/(std::complex<real_t>(1.,0.) + (shifts[index] - shifts.back())* alpha);
#pragma omp parallel for
			for (int site = 0; site < r.localsize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solutions[index][site][mu] = solutions[index][site][mu] + f[index]*alpha*r[site][mu];
				}
			}
			//solutions[index].updateHalo();//not needed anymore
		}

		//Set the last solution
#pragma omp parallel for
		for (int site = 0; site < r.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				r[site][mu] = r[site][mu] - alpha*p[site][mu];
			}
		}
		//r.updateHalo(); not needed anymore

		long_real_t error = AlgebraUtils::squaredNorm(r);
		if (error < epsilon) {
			for (unsigned int index = 0; index < shifts.size(); ++index) {
				original_solutions[index] = solutions[index];
			}
#ifdef TEST_PAPI_SPEED
			if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
			if (isOutputProcess()) std::cout << "MFLOPS for Single Krilov Space Solver: " << mflops << " MFLOPS. " << std::endl;
#endif
			if (isOutputProcess()) std::cout << "MMMRMultishiftSolver::Convergence in " << i << " steps" << std::endl;
			delete[] f;
			return flag;
		}
	}

	for (unsigned int index = 0; index < shifts.size(); ++index) {
		original_solutions[index] = solutions[index];
	}

	long_real_t last_error = AlgebraUtils::squaredNorm(r);

	if (isOutputProcess()) std::cout << "Failure in finding convergence in MMMRMultishiftSolver! " << dirac->getKappa() << " " << last_error << std::endl;
	delete[] f;
	return false;
}

} /* namespace Update */
