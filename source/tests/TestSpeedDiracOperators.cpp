#include "TestSpeedDiracOperators.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_operators/SquareDiracWilsonOperator.h"
#include "dirac_operators/SquareImprovedDiracWilsonOperator.h"
#include "dirac_operators/DiracWilsonOperator.h"
#include "dirac_operators/ImprovedDiracWilsonOperator.h"
#include "dirac_operators/BasicDiracWilsonOperator.h"
#include "dirac_operators/SquareBlockDiracWilsonOperator.h"
#include "dirac_operators/ComplementBlockDiracOperator.h"
#include "dirac_operators/SquareComplementBlockDiracWilsonOperator.h"
#include "dirac_operators/SquareComplementBlockDiracOperator.h"
#include "dirac_operators/SquareTwistedDiracOperator.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "utils/ToString.h"
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

TestSpeedDiracOperators::TestSpeedDiracOperators() : LatticeSweep() { }

TestSpeedDiracOperators::~TestSpeedDiracOperators() { }

void TestSpeedDiracOperators::execute(environment_t& environment) {
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();

	int numberTests;
	try {
		numberTests = environment.configurations.get<unsigned int>("TestSpeedDiracOperators::number_multiplication_test_speed");
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
	SquareImprovedDiracWilsonOperator* squareImprovedDiracWilsonOperator = new SquareImprovedDiracWilsonOperator();
	squareImprovedDiracWilsonOperator->setKappa(0.1);
	squareImprovedDiracWilsonOperator->setCSW(1.);
	squareImprovedDiracWilsonOperator->setLattice(environment.getFermionLattice());
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
		squareImprovedDiracWilsonOperator->multiply(test2,test1);
	}
	clock_gettime(CLOCK_REALTIME, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	if (isOutputProcess()) std::cout << "Timing for SquareImprovedDiracWilsonOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
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
	delete diracWilsonOperator;
	delete improvedDiracWilsonOperator;
	delete basicDiracWilsonOperator;
	delete squareImprovedDiracWilsonOperator;
#ifdef TEST_PAPI_SPEED
	PAPI_shutdown();
#endif
	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
}

void TestSpeedDiracOperators::registerParameters(po::options_description& desc) {
	desc.add_options()
		("TestSpeedDiracOperators::number_multiplication_test_speed", po::value<unsigned int>(), "How many multiplications should I use in the tests?")
		;
}

} /* namespace Update */
