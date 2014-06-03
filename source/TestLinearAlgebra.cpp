/*
 * TestLinearAlgebra.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: spiem_01
 */

#include "TestLinearAlgebra.h"
#include "BiConjugateGradient.h"
#include "AlgebraUtils.h"
#include "dirac_operators/SquareDiracWilsonOperator.h"
#include "dirac_operators/SquareImprovedDiracWilsonOperator.h"
#include "dirac_operators/DiracWilsonOperator.h"
#include "dirac_operators/ImprovedDiracWilsonOperator.h"
#include "dirac_operators/BasicDiracWilsonOperator.h"
#include "dirac_operators/SquareBlockDiracWilsonOperator.h"
#include "dirac_operators/SquareComplementBlockDiracWilsonOperator.h"
#include "dirac_operators/SquareTwistedDiracOperator.h"
#include "DeflationInverter.h"
#include <vector>
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

double diffclock(clock_t clock1,clock_t clock2) {
	double diffticks=clock2-clock1;
	double diffms=(diffticks*1000.)/CLOCKS_PER_SEC;
	return diffms;
}

namespace Update {

TestLinearAlgebra::TestLinearAlgebra() : LatticeSweep() { }

TestLinearAlgebra::~TestLinearAlgebra() { }

void TestLinearAlgebra::execute(environment_t& environment) {
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();

	/*{
		extended_dirac_vector_t source, tmp1, tmp2, tmp3, tmp4;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);

		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		DiracOperator* squareDiracWilsonOperator = new SquareDiracWilsonOperator();
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());
		squareDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		biConjugateGradient->solve(squareDiracWilsonOperator, source, tmp1);
		if (isOutputProcess())	std::cout << "With the square of the dirac wilson operator: " << biConjugateGradient->getLastSteps() << std::endl;
		DiracWilsonOperator* diracWilsonOperator = new DiracWilsonOperator();
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		diracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		biConjugateGradient->solve(diracWilsonOperator, source, tmp2);
		unsigned int steps = biConjugateGradient->getLastSteps();
		biConjugateGradient->solve(diracWilsonOperator, tmp2, tmp3);
		long_real_t diffnorm = AlgebraUtils::differenceNorm(tmp1, tmp3);
		if (isOutputProcess()) {
			std::cout << "With the normal dirac wilson operator: " << steps + biConjugateGradient->getLastSteps() << std::endl;
			std::cout << "Square inverse vs double inverse test: " << diffnorm << std::endl;
		}
		squareDiracWilsonOperator->multiply(tmp4,tmp3);
		diffnorm = AlgebraUtils::differenceNorm(source, tmp4);
		if (isOutputProcess()) {
			std::cout << "Square inverse test : " << diffnorm << std::endl;
		}


		DiracOperator* squareImprovedDiracWilsonOperator = new SquareImprovedDiracWilsonOperator(environment.getFermionLattice(), environment.configurations.get<double>("kappa"), environment.configurations.get<double>("csw"));
		biConjugateGradient->solve(squareImprovedDiracWilsonOperator, source, tmp1);
		if (isOutputProcess())	std::cout << "With the square of the improved dirac wilson operator: " << biConjugateGradient->getLastSteps() << std::endl;
		ImprovedDiracWilsonOperator* improvedDiracWilsonOperator = new ImprovedDiracWilsonOperator();
		improvedDiracWilsonOperator->setLattice(environment.getFermionLattice());
		improvedDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		improvedDiracWilsonOperator->setCSW(environment.configurations.get<double>("csw"));
		biConjugateGradient->solve(improvedDiracWilsonOperator, source, tmp2);
		steps = biConjugateGradient->getLastSteps();
		biConjugateGradient->solve(improvedDiracWilsonOperator, tmp2, tmp3);
		diffnorm = AlgebraUtils::differenceNorm(tmp1, tmp3);
		if (isOutputProcess()) {
			std::cout << "With the improved dirac wilson operator: " << steps + biConjugateGradient->getLastSteps() << std::endl;
			std::cout << "Square inverse vs double inverse test: " << diffnorm << std::endl;
		}
		squareImprovedDiracWilsonOperator->multiply(tmp4,tmp3);
		diffnorm = AlgebraUtils::differenceNorm(source, tmp4);
		if (isOutputProcess()) {
			std::cout << "Square inverse test : " << diffnorm << std::endl;
		}

		delete diracWilsonOperator;
		delete squareDiracWilsonOperator;
		delete improvedDiracWilsonOperator;
		delete squareImprovedDiracWilsonOperator;
		delete biConjugateGradient;
	}*/

	//Gamma5 test
	{
		reduced_dirac_vector_t test1, test2, test3;
		AlgebraUtils::generateRandomVector(test1);
		DiracWilsonOperator* disastro = new DiracWilsonOperator();
		disastro->setLattice(environment.getFermionLattice());
		disastro->setKappa(0.15);
		disastro->multiply(test3, test1);
		disastro->setGamma5(false);
		disastro->multiply(test2, test1);
		long_real_t g5test = norm(AlgebraUtils::dot(test3,test3) - AlgebraUtils::gamma5dot(test3,test2));
		if (isOutputProcess()) std::cout << "Gamma5 test on DiracWilsonOperator: " << g5test << std::endl;
		ImprovedDiracWilsonOperator* disastro2 = new ImprovedDiracWilsonOperator();
		disastro2->setLattice(environment.getFermionLattice());
		disastro2->setKappa(0.15);
		disastro2->setCSW(1.);
		disastro2->multiply(test3, test1);
		disastro2->setGamma5(false);
		disastro2->multiply(test2, test1);
		g5test = norm(AlgebraUtils::dot(test3,test3) - AlgebraUtils::gamma5dot(test3,test2));
		if (isOutputProcess()) std::cout << "Gamma5 test on ImprovedDiracWilsonOperator: " << g5test << std::endl;
	}

	/*{
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.00000000001);
		reduced_dirac_vector_t source, tmp1, tmp2, tmp3, tmp4, tmp5;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);
		
		DiracOperator* squareDiracWilsonOperator = new SquareDiracWilsonOperator();
		squareDiracWilsonOperator->setKappa(0.20);
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());
		
		biConjugateGradient->solve(squareDiracWilsonOperator, source, tmp1);
		int stepsMax = biConjugateGradient->getLastSteps();
		std::cout << "Undeflated number of steps: " << stepsMax << std::endl;

		DeflationInverter* deflationInverter = new DeflationInverter();
		deflationInverter->setPrecision(0.00000000001);
		deflationInverter->setBasisDimension(10);
		deflationInverter->setBlockDivision(4);
		//deflationInverter->generateBasis(squareDiracWilsonOperator);
		deflationInverter->solve(squareDiracWilsonOperator,source,tmp2);
		double normError1 = AlgebraUtils::differenceNorm(tmp1,tmp2);
		squareDiracWilsonOperator->multiply(tmp3, tmp2);
		double normError2 = AlgebraUtils::differenceNorm(tmp3,source);
		squareDiracWilsonOperator->multiply(tmp4, tmp1);
		double normError3 = AlgebraUtils::differenceNorm(tmp4,source);
		if (isOutputProcess()) std::cout << "Ridiamo di gusto: " << normError1 << " " << normError2 << " " << " " << normError3 << " " << deflationInverter->getLastSteps() << std::endl;
		if (isOutputProcess()) std::cout << "Zum beispiel: " << tmp2[5][3] << std::endl << tmp1[5][3] << std::endl;
		delete biConjugateGradient;
		delete squareDiracWilsonOperator;
		delete deflationInverter;
		//delete squareBlockDiracWilsonOperator;
	}*/

	/*{
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.00000000001);
		reduced_dirac_vector_t source, tmp1, tmp2, tmp3, tmp4, tmp5;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);
		DiracOperator* squareDiracWilsonOperator = new SquareDiracWilsonOperator();
		squareDiracWilsonOperator->setKappa(0.205);
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());

		biConjugateGradient->solve(squareDiracWilsonOperator, source, tmp1);
		int stepsMax = biConjugateGradient->getLastSteps();
		std::cout << "Undeflated number of steps: " << stepsMax << std::endl;

		DeflationInverter* deflationInverter = new DeflationInverter();
		deflationInverter->setPrecision(0.00000000001);
		deflationInverter->setBasisDimension(50);
		deflationInverter->setBlockDivision(4);
		//deflationInverter->generateBasis(squareDiracWilsonOperator);
		//deflationInverter->solve(squareDiracWilsonOperator,source,tmp2);
		double normError1 = AlgebraUtils::differenceNorm(tmp1,tmp2);
		squareDiracWilsonOperator->multiply(tmp3, tmp2);
		double normError2 = AlgebraUtils::differenceNorm(tmp3,source);
		squareDiracWilsonOperator->multiply(tmp4, tmp1);
		double normError3 = AlgebraUtils::differenceNorm(tmp4,source);
		if (isOutputProcess()) std::cout << "Ridiamo di gusto: " << normError1 << " " << normError2 << " " << " " << normError3 << " " << deflationInverter->getLastSteps() << std::endl;
		if (isOutputProcess()) std::cout << "Zum beispiel: " << tmp2[5][3] << std::endl << tmp1[5][3] << std::endl;
		delete biConjugateGradient;
		delete squareDiracWilsonOperator;
		//delete deflationInverter;
	}*/

	//Determinant test
	{
		reduced_dirac_vector_t temp1, temp2, result;
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.00000000001);

		SquareImprovedDiracWilsonOperator* squareImprovedDiracWilsonOperator = new SquareImprovedDiracWilsonOperator();
		squareImprovedDiracWilsonOperator->setCSW(1.);
		squareImprovedDiracWilsonOperator->setKappa(0.18);
		squareImprovedDiracWilsonOperator->setLattice(environment.getFermionLattice());

		SquareTwistedDiracOperator* squareTwistedDiracOperator = new SquareTwistedDiracOperator();
		squareTwistedDiracOperator->setDiracOperator(squareImprovedDiracWilsonOperator);
		squareTwistedDiracOperator->setTwist(0.0005);
		squareTwistedDiracOperator->setKappa(0.18);
		squareTwistedDiracOperator->setLattice(environment.getFermionLattice());

		std::vector<long_real_t> stochastic_estimates;

		for (int i = 0; i < 30; ++i) {
			AlgebraUtils::generateRandomGaussianVector(temp1);
			biConjugateGradient->solve(squareTwistedDiracOperator, temp1, temp2);
			squareImprovedDiracWilsonOperator->multiply(result, temp2);
		
			std::complex<long_real_t> saturno = exp(AlgebraUtils::dot(result, temp1) - AlgebraUtils::dot(result,result));
			stochastic_estimates.push_back(real(saturno));

			if (isOutputProcess()) std::cout << "Stima stocastica del determinante: " << saturno << std::endl;
		}

		long_real_t average = 0.;
		long_real_t sd = 0.;
		for (std::vector<long_real_t>::iterator i = stochastic_estimates.begin(); i != stochastic_estimates.end(); ++i) {
			average += *i;
			sd += (*i)*(*i);
		}
		sd = sd/stochastic_estimates.size();
		average = average/stochastic_estimates.size();

		if (isOutputProcess()) {
			std::cout << "Average stochastic estimation of the reweighting factor: " << pow(average,0.25) << std::endl;
		}
	}

	//Hermitian test
	{
		reduced_dirac_vector_t test1, test2, test3, test4;
		AlgebraUtils::generateRandomVector(test1);
		AlgebraUtils::generateRandomVector(test2);
		DiracWilsonOperator* disastro = new DiracWilsonOperator();
		disastro->setLattice(environment.getFermionLattice());
		disastro->setKappa(0.1);
		disastro->multiply(test3, test1);
		disastro->multiply(test4, test2);
		long_real_t htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "Hermitian test on DiracWilsonOperator: " << htest << std::endl;
		AlgebraUtils::generateRandomVector(test1);
		AlgebraUtils::generateRandomVector(test2);
		ImprovedDiracWilsonOperator* disastro2 = new ImprovedDiracWilsonOperator();
		disastro2->setLattice(environment.getFermionLattice());
		disastro2->setKappa(0.1);
		disastro2->setCSW(1.);
		disastro2->multiply(test3,test1);
		disastro2->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "Hermitian test on ImprovedDiracWilsonOperator: " << htest << std::endl;
		BlockDiracWilsonOperator* disastro3 = new BlockDiracWilsonOperator();
		disastro3->setLattice(environment.getFermionLattice());
		disastro3->setKappa(0.1);
		disastro3->multiply(test3,test1);
		disastro3->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "Hermitian test on BlockDiracWilsonOperator: " << htest << std::endl;
		SquareBlockDiracWilsonOperator* disastro4 = new SquareBlockDiracWilsonOperator();
		disastro4->setLattice(environment.getFermionLattice());
		disastro4->setKappa(0.1);
		disastro4->multiply(test3,test1);
		disastro4->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "Hermitian test on SquareBlockDiracWilsonOperator: " << htest << std::endl;
		SquareComplementBlockDiracWilsonOperator* disastro5 = new SquareComplementBlockDiracWilsonOperator();
		disastro5->setLattice(environment.getFermionLattice());
		disastro5->setKappa(0.1);
		disastro5->multiply(test3,test1);
		disastro5->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "Hermitian test on SquareComplementBlockDiracWilsonOperator: " << htest << std::endl;
		delete disastro;
		delete disastro2;
		delete disastro3;
		delete disastro4;
		delete disastro5;
	}

	//Zero test csw
	{
		reduced_dirac_vector_t test1, test2, test3;
		AlgebraUtils::generateRandomVector(test1);
		DiracWilsonOperator* diracWilsonOperator = new DiracWilsonOperator();
		diracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		ImprovedDiracWilsonOperator* improvedDiracWilsonOperator = new ImprovedDiracWilsonOperator();
		improvedDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		improvedDiracWilsonOperator->setCSW(0.0);
		improvedDiracWilsonOperator->setLattice(environment.getFermionLattice());
		diracWilsonOperator->multiplyAdd(test2,test1,test1,0.0);
		improvedDiracWilsonOperator->multiplyAdd(test3,test1,test1,0.0);
		long_real_t zerotest = AlgebraUtils::differenceNorm(test2,test3);
		if (isOutputProcess()) std::cout << "Zero csw test on DiracWilsonOperator (add): " << zerotest << std::endl;
		diracWilsonOperator->multiply(test2,test1);
		improvedDiracWilsonOperator->multiply(test3,test1);
		zerotest = AlgebraUtils::differenceNorm(test2,test3);
		if (isOutputProcess()) std::cout << "Zero csw test on DiracWilsonOperator (no add): " << zerotest << std::endl;
		improvedDiracWilsonOperator->multiply(test2,test1);
		improvedDiracWilsonOperator->multiplyAdd(test3,test1,test1,0.0);
		zerotest = AlgebraUtils::differenceNorm(test2,test3);
		if (isOutputProcess()) std::cout << "Zero alpha test on ImprovedWilsonOperator: " << zerotest << std::endl;
		BasicDiracWilsonOperator* basicDiracWilsonOperator = new BasicDiracWilsonOperator();
		basicDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		basicDiracWilsonOperator->setLattice(environment.getFermionLattice());
		basicDiracWilsonOperator->multiply(test3,test1);
		zerotest = AlgebraUtils::differenceNorm(test2,test3);
		if (isOutputProcess()) std::cout << "Zero test on DiracWilsonOperator vs BasicDiracWilsonOperator (no add): " << zerotest << std::endl;
		delete diracWilsonOperator;
		delete improvedDiracWilsonOperator;
		delete basicDiracWilsonOperator;
	}

	{
		int numberTests;
		try {
			numberTests = environment.configurations.get<unsigned int>("number_multiplication_test_speed");
		} catch (NotFoundOption& e) {
			numberTests = 300;
		}
		DiracWilsonOperator* diracWilsonOperator = new DiracWilsonOperator();
		diracWilsonOperator->setKappa(0.1);
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		BasicDiracWilsonOperator* basicDiracWilsonOperator = new BasicDiracWilsonOperator();
		basicDiracWilsonOperator->setKappa(0.1);
		basicDiracWilsonOperator->setLattice(environment.getFermionLattice());
		ImprovedDiracWilsonOperator* improvedDiracWilsonOperator = new ImprovedDiracWilsonOperator();
		improvedDiracWilsonOperator->setKappa(0.1);
		improvedDiracWilsonOperator->setCSW(1.);
		improvedDiracWilsonOperator->setLattice(environment.getFermionLattice());
		reduced_dirac_vector_t test1, test2, test3, test4, test5, test6, test7, test8;
		AlgebraUtils::generateRandomVector(test1);
		AlgebraUtils::generateRandomVector(test3);
		AlgebraUtils::generateRandomVector(test5);
		AlgebraUtils::generateRandomVector(test7);
#ifdef TEST_PAPI_SPEED
		float real_time, proc_time, mflops;
		long long flpins;
		int retval;
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
#endif
		struct timespec start, finish;
		double elapsed;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			diracWilsonOperator->multiply(test2,test1);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
		if (isOutputProcess()) std::cout << "MFLOPS for DiracWilsonOperator: " << mflops << " MFLOPS. " << std::endl;
#endif
		if (isOutputProcess()) std::cout << "Timing for DiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			diracWilsonOperator->multiply(test2,test4,test1,test3);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "Timing for double DiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			diracWilsonOperator->multiply(test2,test4,test6,test8,test1,test3,test5,test7);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "Timing for quad DiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			basicDiracWilsonOperator->multiply(test2,test1);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "Timing for BasicDiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
#endif
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			improvedDiracWilsonOperator->multiply(test2,test1);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
		if (isOutputProcess()) std::cout << "MFLOPS for ImprovedDiracWilsonOperator: " << mflops << " MFLOPS. " << std::endl;
#endif
		if (isOutputProcess()) std::cout << "Timing for ImprovedDiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms." << std::endl;
		
		std::complex<long_real_t> result_fast;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
#endif
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			result_fast = AlgebraUtils::dot(test3,test1);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
		if (isOutputProcess()) std::cout << "MFLOPS for AlgebraUtils::dot: " << mflops << " MFLOPS. " << std::endl;
#endif
		if (isOutputProcess()) std::cout << "Timing for AlgebraUtils::dot: " << (elapsed*1000)/numberTests << " ms." << std::endl;
		
		/*std::complex<long_real_t> result_slow;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			result_slow = AlgebraUtils::slow_dot(test3,test1);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "Timing for AlgebraUtils::slow_dot: " << (elapsed*1000)/numberTests << " ms." << std::endl;
		if (isOutputProcess()) std::cout << result_fast - result_slow << std::endl;*/
		
		delete diracWilsonOperator;
		delete improvedDiracWilsonOperator;
		delete basicDiracWilsonOperator;
	}
#ifdef TEST_PAPI_SPEED
	PAPI_shutdown();
#endif
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
}

} /* namespace Update */
