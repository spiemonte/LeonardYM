/*
 * TestLinearAlgebra.cpp
 *
 *  Created on: Jul 3, 2012
 *      Author: spiem_01
 */

#include "TestLinearAlgebra.h"
#include "inverters/BiConjugateGradient.h"
#include "inverters/GMRESR.h"
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
#include "dirac_operators/MultiGridOperator.h"
#include "inverters/DeflationInverter.h"
#include "dirac_functions/Polynomial.h"
#include "utils/ToString.h"
#include <vector>


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
		biConjugateGradient->setMaximumSteps(8000);
		DiracOperator* squareDiracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareDiracOperator->setLattice(environment.getFermionLattice());
		
		biConjugateGradient->solve(squareDiracOperator, source, tmp1);
		if (isOutputProcess())	std::cout << "With the square of the dirac wilson operator: " << biConjugateGradient->getLastSteps() << std::endl;
		
		DiracOperator* diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
		diracOperator->setLattice(environment.getFermionLattice());

		diracOperator->setGamma5(false);
		AlgebraUtils::gamma5(source);
		biConjugateGradient->solve(diracOperator, source, tmp2);
		
		
		unsigned int steps = biConjugateGradient->getLastSteps();
		AlgebraUtils::gamma5(tmp2);
		biConjugateGradient->solve(diracOperator, tmp2, tmp3);
		
		long_real_t diffnorm = AlgebraUtils::differenceNorm(tmp1, tmp3);
		if (isOutputProcess()) {
			std::cout << "With the normal dirac wilson operator: " << steps + biConjugateGradient->getLastSteps() << std::endl;
			std::cout << "Square inverse vs double inverse test: " << diffnorm << std::endl;
		}
		squareDiracOperator->multiply(tmp4,tmp3);
		AlgebraUtils::gamma5(source);
		diffnorm = AlgebraUtils::differenceNorm(source, tmp4);
		if (isOutputProcess()) {
			std::cout << "Square inverse test : " << diffnorm << std::endl;
		}

		delete diracOperator;
		delete squareDiracOperator;
delete biConjugateGradient;
}*/


	{
		/*DiracOperator* squareDiracWilsonOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());*/

		DiracOperator* diracWilsonOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		diracWilsonOperator->setGamma5(false);

		std::vector<unsigned int> blockSize(4);
		blockSize[0] = 4;
		blockSize[1] = 4;
		blockSize[2] = 4;
		blockSize[3] = 4;

		BlockDiracOperator* blackBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Black);
		blackBlockDiracOperator->setLattice(environment.getFermionLattice());
		blackBlockDiracOperator->setGamma5(false);
		blackBlockDiracOperator->setBlockSize(blockSize);
		
		BlockDiracOperator* redBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Red);
		redBlockDiracOperator->setLattice(environment.getFermionLattice());
		redBlockDiracOperator->setGamma5(false);
		redBlockDiracOperator->setBlockSize(blockSize);

		TwistedDiracOperator* twistedDiracOperator = new TwistedDiracOperator();
		twistedDiracOperator->setDiracOperator(diracWilsonOperator);
		twistedDiracOperator->setLattice(environment.getFermionLattice());
		//twistedDiracOperator->setTwist(0.0001);
		
		ComplementBlockDiracOperator* K = new ComplementBlockDiracOperator(diracWilsonOperator, redBlockDiracOperator, blackBlockDiracOperator);
		K->setMaximumSteps(205);
		K->setBlockSize(blockSize);
		
		


		
		SAPPreconditioner *preconditioner = new SAPPreconditioner(twistedDiracOperator,K);
		preconditioner->setSteps(7);
		preconditioner->setPrecision(0.00001);

		MultiGridOperator* multiGridOperator = new MultiGridOperator();
		MultiGridProjector* multiGridProjector = new MultiGridProjector();
		multiGridOperator->setDiracOperator(diracWilsonOperator);
		
		int basisDimension = 20;
		MultiGridVectorLayout::setBasisDimension(basisDimension);
		MultiGridVectorLayout::tBlockSize = 4;
		MultiGridVectorLayout::initialize();

		BlockBasis test(basisDimension);
		reduced_dirac_vector_t zeroVector;
		AlgebraUtils::setToZero(zeroVector);
		reduced_dirac_vector_t randomVector;
		GMRESR* gmres_inverter = new GMRESR();
		gmres_inverter->setPrecision(0.0000000001);
		gmres_inverter->setMaximumSteps(1);

		struct timespec start, finish;
		double elapsed;
		clock_gettime(CLOCK_REALTIME, &start);
		
		for (int i = 0; i < basisDimension; i += 1) {
			//We start with a random vector
			AlgebraUtils::generateRandomVector(randomVector);
			AlgebraUtils::normalize(randomVector);

			//We give random vector as initial guess to solve the omogeneous system
			gmres_inverter->solve(diracWilsonOperator, zeroVector, test[i], preconditioner, &randomVector);
			//test[i+1] = test[i];
			//AlgebraUtils::gamma5(test[i+1]);
			/*diracWilsonOperator->multiply(randomVector,test[i]);
			if (isOutputProcess()) std::cout << "MultiGrid::Deficit for the vector " << i << ": " << AlgebraUtils::squaredNorm(randomVector)/AlgebraUtils::squaredNorm(test[i]) << "" << std::endl;*/
		}

		test.orthogonalize();

		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Multigrid basis constructed in: " << (elapsed) << " s."<< std::endl;

		for (int i = 0; i < basisDimension; ++i) {
			multiGridOperator->addVector(test[i]);
			multiGridProjector->addVector(test[i]);
		}

		int numberTests;
		try {
			numberTests = environment.configurations.get<unsigned int>("number_multiplication_test_speed");
		} catch (NotFoundOption& e) {
			numberTests = 300;
		}

		multigrid_vector_t tmp1, tmp2;
		for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) tmp2[i] = std::complex<real_t>(0.2-i*i,0.4+i);

		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			multiGridOperator->multiply(tmp1,tmp2);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
		if (isOutputProcess()) std::cout << "MFLOPS for MultiGridOperator: " << mflops << " MFLOPS. " << std::endl;
#endif
		if (isOutputProcess()) std::cout << "Timing for MultiGridOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;

		
		MultiGridBiConjugateGradientSolver* biMgSolver = new MultiGridBiConjugateGradientSolver();
		biMgSolver->setMaximumSteps(35);
		biMgSolver->setPrecision(0.00000000001);

		multigrid_vector_t solution_hat, source_hat;


		int conjugateSpaceDimension = 15;
		reduced_dirac_vector_t source, r, solution, c[conjugateSpaceDimension], u[conjugateSpaceDimension], mg_inverse, source_sap;

		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);

		r = source;
		AlgebraUtils::setToZero(solution);

		clock_gettime(CLOCK_REALTIME, &start);

		for (int k = 0; k < 535; ++k) {
			{
				multiGridOperator->setDiracOperator(diracWilsonOperator);

				multiGridProjector->apply(source_hat,r);
				biMgSolver->solve(multiGridOperator, source_hat, solution_hat);

				//multiGridOperator->multiply(source_test,solution_hat);

				multiGridProjector->apply(mg_inverse,solution_hat);

				diracWilsonOperator->multiply(source_sap,mg_inverse);
#pragma omp parallel for
				for (int site = 0; site < solution.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						source_sap[site][mu] = r[site][mu] - source_sap[site][mu];
					}
				}
				
				preconditioner->multiply(u[(k % conjugateSpaceDimension)],source_sap);
#pragma omp parallel for
				for (int site = 0; site < solution.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu] + mg_inverse[site][mu];
					}
				}

				//preconditioner->multiply(u[(k % conjugateSpaceDimension)],r);
				
				//preconditioner->multiply(tmp1,r);

				//diracWilsonOperator->multiply(tmp2,tmp1);
				//diracWilsonOperator->multiply(tmp2,u[(k % conjugateSpaceDimension)]);
				//std::cout << "Giusto per: " << AlgebraUtils::differenceNorm(r,tmp2) << " " << AlgebraUtils::differenceNorm(r,tmp3) << std::endl;
			}
			//u[(k % conjugateSpaceDimension)] = r;

			diracWilsonOperator->multiply(c[(k % conjugateSpaceDimension)],u[(k % conjugateSpaceDimension)]);
			
			for (int i = 0; i < (k % conjugateSpaceDimension); ++i) {
				std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(c[i],c[(k % conjugateSpaceDimension)]));
#pragma omp parallel for
				for (int site = 0; site < solution.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						c[(k % conjugateSpaceDimension)][site][mu] = c[(k % conjugateSpaceDimension)][site][mu] - alpha*c[i][site][mu];
						u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu] - alpha*u[i][site][mu];
					}
				}
				
			}
			
			real_t norm = sqrt(AlgebraUtils::squaredNorm(c[(k % conjugateSpaceDimension)]));
#pragma omp parallel for
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					c[(k % conjugateSpaceDimension)][site][mu] = c[(k % conjugateSpaceDimension)][site][mu]/norm;
					u[(k % conjugateSpaceDimension)][site][mu] = u[(k % conjugateSpaceDimension)][site][mu]/norm;
				}
			}
			
			std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(c[(k % conjugateSpaceDimension)],r));
#pragma omp parallel for
			for (int site = 0; site < solution.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = solution[site][mu] + alpha*u[(k % conjugateSpaceDimension)][site][mu];
					r[site][mu] = r[site][mu] - alpha*c[(k % conjugateSpaceDimension)][site][mu];
				}
			}

			long_real_t error = AlgebraUtils::squaredNorm(r);
			if (error < 0.00000000001) {
				//stepsMax = j;
				break;
			}
			else if (isOutputProcess()) std::cout << "Residual norm at step " << k << ": " << error << std::endl;

		}

		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Multigrid inversion done in: " << (elapsed) << " s."<< std::endl;

		gmres_inverter->setMaximumSteps(10000);
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.00000000001);
		biConjugateGradient->setMaximumSteps(100000);

		clock_gettime(CLOCK_REALTIME, &start);
		biConjugateGradient->solve(diracWilsonOperator, source, solution);
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		
		if (isOutputProcess()) std::cout << "TestLinearAlgebra:: vs standard inversion done in: " << (elapsed) << " s."<< std::endl;
		
	}
#ifdef TESTALL
	{
		DiracOperator* squareDiracWilsonOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());

		DiracOperator* diracWilsonOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		diracWilsonOperator->setGamma5(false);

		MultiGridOperator* multiGridOperator = new MultiGridOperator();
		MultiGridProjector* multiGridProjector = new MultiGridProjector();
		multiGridOperator->setDiracOperator(squareDiracWilsonOperator);

		MultiGridVectorLayout::initialize();

		BlockBasis test(70);
		reduced_dirac_vector_t zeroVector;
		reduced_dirac_vector_t randomVector;

		ConjugateGradient* conjugateGradient = new ConjugateGradient();

		std::vector< complex > roots;
		roots.push_back(std::complex<real_t>(0.81480744985345959026,-0.73998964487862817075));
		roots.push_back(std::complex<real_t>(0.81480744985345959026,0.73998964487862817075));
		roots.push_back(std::complex<real_t>(2.7411901619266698148,-1.1682614363852250770));
		roots.push_back(std::complex<real_t>(2.7411901619266698148,1.1682614363852250770));
		roots.push_back(std::complex<real_t>(4.937010504853446515,-0.967517145569813877));
		roots.push_back(std::complex<real_t>(4.937010504853446515,0.967517145569813877));
		roots.push_back(std::complex<real_t>(6.51157966762550890,-0.29309099771391338));
		roots.push_back(std::complex<real_t>(6.51157966762550890,0.29309099771391338));
		
		Polynomial polynomialPreconditioner;
		polynomialPreconditioner.setRoots(roots);
		polynomialPreconditioner.setScaling(0.4046234347673221166251352);

		for (int i = 0; i < 70; i += 2) {
			//We start with a random vector
			AlgebraUtils::generateRandomVector(randomVector);

			conjugateGradient->setMaximumSteps(290);
			conjugateGradient->setPrecision(0.00000000001);
			//We give random vector as initial guess to solve the omogeneous system
			conjugateGradient->solve(squareDiracWilsonOperator,zeroVector,test[i],&randomVector);
			/*for (int k = 0; k < 4; ++k) {
				randomVector = test[i];
				AlgebraUtils::normalize(randomVector);
				polynomialPreconditioner.evaluate(squareDiracWilsonOperator,test[i],randomVector);
			}*/
			//test[i] = randomVector;
			test[i+1] = test[i];
			AlgebraUtils::gamma5(test[i+1]);
			squareDiracWilsonOperator->multiply(randomVector,test[i]);
			std::cout << "MultiGrid::Deficit for the vector " << i << ": " << AlgebraUtils::squaredNorm(randomVector)/AlgebraUtils::squaredNorm(test[i]) << "" << std::endl;
		}

		test.orthogonalize();

		for (int i = 0; i < 70; ++i) {
			multiGridOperator->addVector(test[i]);
			multiGridProjector->addVector(test[i]);
		}

		//std::cout << toString(multiGridOperator->asMatrix()) << std::endl;



		/*{
			multigrid_vector_t test1, test2, test3, test4;
			for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) {
				test1[i] = std::complex<real_t>(0.2+i*i,0.4-i);
				test2[i] = std::complex<real_t>(0.2-i*i,0.4+i);
			}
			
			multiGridOperator->multiply(test3,test1);
			multiGridOperator->multiply(test4,test2);

			long_real_t norm1 = 0., norm2 = 0.;
#pragma omp parallel for reduction(+:norm1,norm2)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				norm1 += real(conj(test2[m])*test3[m]);
				norm2 += real(conj(test1[m])*test4[m]);
			}

			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on SquareBlockDiracWilsonOperator: " << norm1 - norm2 << std::endl;
		}

		{
			multigrid_vector_t test1, test2;
			reduced_dirac_vector_t tmp3;
			AlgebraUtils::generateRandomVector(tmp3);
			multiGridProjector->apply(test1,tmp3);
			multiGridProjector->apply(tmp3,test1);
			multiGridProjector->apply(test2,tmp3);

			long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				norm += real(conj(test2[m] - test1[m])*(test2[m] - test1[m]));
			}

			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Double projection test: " << norm << std::endl;
		}

		{
			multigrid_vector_t test1, test2;
			reduced_dirac_vector_t tmp3;
			AlgebraUtils::generateRandomVector(tmp3);
			multiGridProjector->apply(test1,tmp3);
			multiGridOperator->multiply(test2,test1);

			long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
			for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
				norm += real(conj(test2[m] - test1[m])*(test2[m] - test1[m]));
			}

			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Identity test: " << norm << std::endl;
		}*/



		int stepsMax = 15;

		reduced_dirac_vector_t source, guess, tmp, tmp3, tmp4, tmp5, tmp6, tmp7, Kb, p, chi, r, r_hat, solution, update;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);
		solution = source;

		conjugateGradient->setPrecision(0.0000000000001);
		conjugateGradient->setMaximumSteps(13000);

		struct timespec start, finish;
		double elapsed;
		clock_gettime(CLOCK_REALTIME, &start);

		//conjugateGradient->solve(squareDiracWilsonOperator,source,tmp5);


		std::cout << "TestLinearAlgebra::Undeflated solution in " << conjugateGradient->getLastSteps() << " steps." << std::endl;
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "TestLinearAlgebra:: and in: " << (elapsed) << " s."<< std::endl;

		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.0000000000001);
		biConjugateGradient->setMaximumSteps(100);

		MultiGridConjugateGradientSolver* mgSolver = new MultiGridConjugateGradientSolver();
		mgSolver->setMaximumSteps(100);
		mgSolver->setPrecision(0.0000000000001);

		MultiGridBiConjugateGradientSolver* biMgSolver = new MultiGridBiConjugateGradientSolver();
		biMgSolver->setMaximumSteps(10);
		biMgSolver->setPrecision(0.00000000001);

		solution = source;


		/*{
			reduced_dirac_vector_t tmp1,tmp2,tmp3,tmp4;
			squareDiracWilsonOperator->multiply(tmp1,source);
			leftProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp2,tmp1);
			rightProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp3,source);
			squareDiracWilsonOperator->multiply(tmp4,tmp3);
			std::cout << "Primo test: " << AlgebraUtils::differenceNorm(tmp4,tmp2) << std::endl;

			leftProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp1,source);
			leftProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp2,tmp1);
			std::cout << "Secondo test: " << AlgebraUtils::differenceNorm(tmp1,tmp2) << std::endl;

			rightProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp1,source);
			rightProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, tmp2,tmp1);
			std::cout << "Terzo test: " << AlgebraUtils::differenceNorm(tmp1,tmp2) << std::endl;
		}

		reduced_dirac_vector_t projectedSource;
		leftProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, projectedSource,source);

		conjugateGradient->setMaximumSteps(15);
		conjugateGradient->solve(deflatedDirac,projectedSource,chi);*/

		/*r = projectedSource;
		AlgebraUtils::setToZero(chi);
		int conjugateSpaceDimension = 15;
		std::complex<real_t> rho[conjugateSpaceDimension];
		for (int j = 0; j < 35; ++j) {
			//Inner conjugate gradient solver

			w[(j % conjugateSpaceDimension)] = r;

			d[(j % conjugateSpaceDimension)] = w[(j % conjugateSpaceDimension)];

			int js = 0;//((j % conjugateSpaceDimension) == 0) ? 0 : 1;
			for (int k = js; k < (j % conjugateSpaceDimension); ++k) {
				std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(w[(j % conjugateSpaceDimension)],ad[k]))/rho[k];
#pragma omp parallel for
				for (int site = 0; site < guess.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						d[(j % conjugateSpaceDimension)][site][mu] = d[(j % conjugateSpaceDimension)][site][mu] - alpha*d[k][site][mu];
					}
				}
			}

			deflatedDirac->multiply(ad[(j % conjugateSpaceDimension)],d[(j % conjugateSpaceDimension)]);
			rho[(j % conjugateSpaceDimension)] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[(j % conjugateSpaceDimension)],ad[(j % conjugateSpaceDimension)]));

			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[(j % conjugateSpaceDimension)],r))/rho[(j % conjugateSpaceDimension)];

#pragma omp parallel for
			for (int site = 0; site < guess.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					chi[site][mu] = chi[site][mu] + proj*d[(j % conjugateSpaceDimension)][site][mu];
					r[site][mu] = r[site][mu] - proj*ad[(j % conjugateSpaceDimension)][site][mu];
				}
			}
			//std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;
			if (AlgebraUtils::squaredNorm(r) < 0.00000000001) {
				stepsMax = j;
				break;
			}
			else std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;

		}*/

		

		/*std::cout << "Deflated convergence (e ridiamo) in " << biConjugateGradient->getLastSteps() << " steps." << std::endl;

		reduced_dirac_vector_t chiProjected;
		rightProjector->apply(squareDiracWilsonOperator, multiGridOperator, multiGridProjector, chiProjected, chi);

		reduced_dirac_vector_t eta;

		multigrid_vector_t solution_hat, source_hat;
		multiGridProjector->apply(source_hat,source);

		mgSolver->solve(multiGridOperator, source_hat, solution_hat);

		multiGridProjector->apply(eta,solution_hat);

#pragma omp parallel for
		for (int site = 0; site < solution.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				solution[site][mu] = eta[site][mu] + chiProjected[site][mu];
			}
		}*/


		extended_fermion_lattice_t blockedLattice = environment.getFermionLattice();
		typedef extended_fermion_lattice_t::Layout Layout;

		int xBlockSize = 8, yBlockSize = 8, zBlockSize = 8, tBlockSize = 8;

#pragma omp parallel for
		for (int site = 0; site < blockedLattice.localsize; ++site) {
			if (Layout::globalIndexX(site) % xBlockSize) set_to_zero(blockedLattice[site][0]);
			if (Layout::globalIndexY(site) % yBlockSize) set_to_zero(blockedLattice[site][1]);
			if (Layout::globalIndexZ(site) % zBlockSize) set_to_zero(blockedLattice[site][2]);
			if (Layout::globalIndexT(site) % tBlockSize) set_to_zero(blockedLattice[site][3]);
		}

		blockedLattice.updateHalo();

		



		BlockDiracOperator* blackBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Black);
		blackBlockDiracOperator->setLattice(environment.getFermionLattice());
		blackBlockDiracOperator->setGamma5(false);
		
		BlockDiracOperator* redBlockDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations, Red);
		redBlockDiracOperator->setLattice(environment.getFermionLattice());
		redBlockDiracOperator->setGamma5(false);

		SquareTwistedDiracOperator* squareTwistedDiracOperator = new SquareTwistedDiracOperator();
		squareTwistedDiracOperator->setDiracOperator(diracWilsonOperator);
		
		TwistedDiracOperator* twistedDiracOperator = new TwistedDiracOperator();
		twistedDiracOperator->setDiracOperator(diracWilsonOperator);

		TwistedDiracOperator* twistedRed = new TwistedDiracOperator();
		twistedRed->setDiracOperator(redBlockDiracOperator);

		TwistedDiracOperator* twistedBlack = new TwistedDiracOperator();
		twistedBlack->setDiracOperator(blackBlockDiracOperator);

		ComplementBlockDiracOperator* K = new ComplementBlockDiracOperator(diracWilsonOperator, redBlockDiracOperator, blackBlockDiracOperator);
		K->setMaximumSteps(15);
		
		SquareComplementBlockDiracOperator* K2 = new SquareComplementBlockDiracOperator(K);
		
		SAPPreconditioner *preconditioner = new SAPPreconditioner(twistedDiracOperator,K);
		preconditioner->setSteps(5);
		preconditioner->setPrecision(0.000000000001);
		

		biConjugateGradient->setMaximumSteps(3000);

		clock_gettime(CLOCK_REALTIME, &start);

		MGRightProjector* rightProjector = new MGRightProjector();
		MGLeftProjector* leftProjector = new MGLeftProjector();
		
		MGDeflatedDirac* deflatedDirac = new MGDeflatedDirac(twistedDiracOperator, leftProjector, multiGridOperator, multiGridProjector);

		reduced_dirac_vector_t projectedSource, z, r_next, z_next, ap;
		multigrid_vector_t solution_hat, source_hat, dsource_hat;


		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);


		int gmresBasis = 3;
		reduced_dirac_vector_t dBasis[gmresBasis], zBasis[gmresBasis], mg_residual, mg_inverse, dmg_inverse, sap_inverse, dsap_inverse, dr_hat;
		r = source;
		AlgebraUtils::setToZero(solution);
		
		for (int j = 0; j < 835; ++j) {
			r_hat = r;
			for (int i = 0; i < gmresBasis; ++i) {
				diracWilsonOperator->setGamma5(false);
				/*multiGridOperator->setDiracOperator(diracWilsonOperator);

				multiGridProjector->apply(source_hat,r_hat);

				biMgSolver->solve(multiGridOperator, source_hat, solution_hat);

				multiGridOperator->multiply(source_test,solution_hat);

				multiGridProjector->apply(mg_inverse,solution_hat);*/

				multiGridOperator->setDiracOperator(diracWilsonOperator);

#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 2; ++mu) {
						tmp5[site][mu] = r_hat[site][mu];
					}
					for (unsigned int mu = 2; mu < 4; ++mu) {
						tmp5[site][mu] = -r_hat[site][mu];
					}
				}

				diracWilsonOperator->multiply(tmp6,tmp5);

#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 2; mu < 4; ++mu) {
						tmp6[site][mu] = -tmp6[site][mu];
					}
				}

				multiGridProjector->apply(source_hat,tmp6);

				//multiGridOperator->multiply(dsource_hat,source_hat);

				multiGridOperator->setDiracOperator(squareDiracWilsonOperator);

				biMgSolver->solve(multiGridOperator, source_hat, solution_hat);

				multiGridProjector->apply(mg_inverse,solution_hat);

				diracWilsonOperator->multiply(dmg_inverse,mg_inverse);
#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						mg_residual[site][mu] = r_hat[site][mu] - dmg_inverse[site][mu];
					}
				}

				K->multiply(sap_inverse, mg_residual);

				diracWilsonOperator->multiply(dsap_inverse,sap_inverse);

				diracWilsonOperator->multiply(dr_hat,r_hat);

				matrix_t B(3,3);
				B(0,0) = AlgebraUtils::dot(dmg_inverse,dmg_inverse);
				B(0,1) = AlgebraUtils::dot(dmg_inverse,dsap_inverse);
				B(0,2) = AlgebraUtils::dot(dmg_inverse,dr_hat);
				B(1,0) = AlgebraUtils::dot(dsap_inverse,dmg_inverse);
				B(1,1) = AlgebraUtils::dot(dsap_inverse,dsap_inverse);
				B(1,2) = AlgebraUtils::dot(dsap_inverse,dr_hat);
				B(2,0) = AlgebraUtils::dot(dr_hat,dmg_inverse);
				B(2,1) = AlgebraUtils::dot(dr_hat,dsap_inverse);
				B(2,2) = AlgebraUtils::dot(dr_hat,dr_hat);

				vector_t V(3);
				V(0) = AlgebraUtils::dot(dmg_inverse,r_hat);
				V(1) = AlgebraUtils::dot(dsap_inverse,r_hat);
				V(2) = AlgebraUtils::dot(dr_hat,r_hat);
#ifndef EIGEN
				vector_t c = inverse(B)*V;
#endif
#ifdef EIGEN
				vector_t c = B.colPivHouseholderQr().solve(V);
#endif				

#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						zBasis[i][site][mu] = c(0)*mg_inverse[site][mu] + c(1)*sap_inverse[site][mu] + c(2)*r_hat[site][mu];
					}
				}

				diracWilsonOperator->multiplyAdd(tmp6,r_hat,r_hat,-1.);
				diracWilsonOperator->multiplyAdd(tmp5,zBasis[i],r_hat,-1.);
				std::cout << "Linear coefficients: {" << c(0) << ", " << c(1) << ", " << c(2) << "}, norm correction: " << AlgebraUtils::squaredNorm(tmp6) << " vs " << AlgebraUtils::squaredNorm(tmp5) << std::endl;

				//K->multiply(zBasis[i],r_hat);
		
				AlgebraUtils::normalize(zBasis[i]);
				diracWilsonOperator->multiply(dBasis[i],zBasis[i]);

				r_hat = dBasis[i];
				
			}

			/*for (int i = 0; i < gmresBasis/2; ++i) {
				diracWilsonOperator->setGamma5(false);

				K->multiply(zBasis[i],r_hat);
				
				diracWilsonOperator->multiplyAdd(tmp6,zBasis[i],r_hat,-1.);
				std::cout << "Giusto per 1: " << AlgebraUtils::squaredNorm(tmp6) << std::endl;
				diracWilsonOperator->multiplyAdd(tmp6,r_hat,r_hat,-1.);
				std::cout << "Giusto per 2: " << AlgebraUtils::squaredNorm(tmp6) << std::endl;
				
				AlgebraUtils::normalize(zBasis[i]);
				diracWilsonOperator->multiply(dBasis[i],zBasis[i]);

				r_hat = dBasis[i];
				
			}

			for (int i = gmresBasis/2; i < gmresBasis; ++i) {
				diracWilsonOperator->setGamma5(false);
				multiGridOperator->setDiracOperator(diracWilsonOperator);

				multiGridProjector->apply(source_hat,dBasis[i-gmresBasis/2]);

				biMgSolver->solve(multiGridOperator, source_hat, solution_hat);

				multiGridOperator->multiply(source_test,solution_hat);

				std::cout << "vediamo chi e' lo stronzo: " << source_test[1] << " - " << source_hat[1] << std::endl;

				multiGridProjector->apply(zBasis[i],solution_hat);
		
				AlgebraUtils::normalize(zBasis[i]);
				diracWilsonOperator->multiply(dBasis[i],zBasis[i]);
			}*/				



			matrix_t littleOperator(gmresBasis,gmresBasis);
			vector_t littleSource(gmresBasis), littleSolution;
			for (int m = 0; m < gmresBasis; ++m) {
				for (int n = 0; n < gmresBasis; ++n) {
					littleOperator(m,n) = AlgebraUtils::dot(dBasis[m],dBasis[n]);
				}
				littleSource(m) = AlgebraUtils::dot(dBasis[m],r);
			}
#ifndef EIGEN
			matrix_t inverseLittleOperator = inverse(littleOperator);
			littleSolution = inverseLittleOperator*littleSource;
#endif
#ifdef EIGEN
			littleSolution = littleOperator.colPivHouseholderQr().solve(littleSource);
#endif

			for (int n = 0; n < gmresBasis; ++n) {
#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						solution[site][mu] -= littleSolution[n]*zBasis[n][site][mu];
					}
				}
			}

			diracWilsonOperator->multiply(tmp5,solution);
#pragma omp parallel for
			for (int site = 0; site < source.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					r[site][mu] = tmp5[site][mu] - source[site][mu];
				}
			}

			std::cout << "GMRES::residual at step: " << j << ": "<< AlgebraUtils::squaredNorm(r) << std::endl;

			//test[(j) % 30] = r;
			//test[(j+1) % 30] = test[(j) % 30];
			//AlgebraUtils::gamma5(test[(j+1) % 30]);
			//test.orthogonalize();
			
		}





		
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;





		
		real_t mu = 0.;

		squareDiracWilsonOperator->multiplyAdd(tmp5,solution,solution,mu*mu);
#pragma omp parallel for
		for (int site = 0; site < source.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				r[site][mu] = source[site][mu] - tmp5[site][mu];
			}
		}

		real_t r_norm = AlgebraUtils::squaredNorm(r);
		std::cout << "TestLinearAlgebra::Residual norm square root inverter " << r_norm << " after " << 0 << " steps." << std::endl;

		if (isOutputProcess()) std::cout << "TestLinearAlgebra:: obtained in: " << (elapsed) << " s."<< std::endl;

		/*int basisIndex = 0;
		for (int j = 0; j < 30; ++j) {

			squareDiracWilsonOperator->multiply(tmp7,solution);
			//Set r to b
#pragma omp parallel for
			for (int site = 0; site < r.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					r[site][mu] = source[site][mu] - tmp7[site][mu];
				}
			}
			
			for (unsigned int i = 0; i < 40; ++i) {
				squareDiracWilsonOperator->multiply(p,r);
				std::complex<real_t> alpha = 0.95*static_cast< std::complex<real_t> >((AlgebraUtils::dot(p,r))/(AlgebraUtils::dot(p,p)));

#pragma omp parallel for
				for (int site = 0; site < r.localsize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						solution[site][mu] = solution[site][mu] + alpha*r[site][mu];
					}
				}

#pragma omp parallel for
				for (int site = 0; site < r.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						r[site][mu] = r[site][mu] - alpha*p[site][mu];
					}
				}

				long_real_t error = AlgebraUtils::squaredNorm(r);
				if (error < 0.00000000001) {
					if (isOutputProcess()) std::cout << "MMMRMultishiftSolver::Convergence in " << i << " steps" << std::endl;
				}
			}
			
			squareDiracWilsonOperator->multiply(tmp5,solution);
#pragma omp parallel for
			for (int site = 0; site < source.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					r[site][mu] = source[site][mu] - tmp5[site][mu];
				}
			}
			
			real_t r_norm = AlgebraUtils::squaredNorm(r);
			real_t norm_next = r_norm;
			std::cout << "Residual norm saturno at step " << j << " " << r_norm << " before MG deflation"<< std::endl;

			for (int i = 0; i < 1; ++i) {
				test[(basisIndex % 7)] = r;
				test.orthogonalize();
				++basisIndex;

				multigrid_vector_t solution_hat, source_hat;
				multiGridProjector->apply(source_hat,r);

				mgSolver->solve(multiGridOperator, source_hat, solution_hat);

				multiGridProjector->apply(tmp6,solution_hat);

				squareDiracWilsonOperator->multiply(tmp5,tmp6);

				std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(tmp5,r)/AlgebraUtils::squaredNorm(tmp5));
				                                                                 
#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						solution[site][mu] = solution[site][mu] + alpha*tmp6[site][mu];
					}
				}


				squareDiracWilsonOperator->multiply(tmp5,solution);
#pragma omp parallel for
				for (int site = 0; site < source.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						r[site][mu] = source[site][mu] - tmp5[site][mu];
					}
				}

				norm_next = AlgebraUtils::squaredNorm(r);
				std::cout << "Residual norm saturno at step " << j << " " << norm_next << " after MG deflation" << std::endl;

				if (norm_next > r_norm && i > 0) {
#pragma omp parallel for
					for (int site = 0; site < source.completesize; ++site) {
						for (unsigned int mu = 0; mu < 4; ++mu) {
							solution[site][mu] = solution[site][mu] - tmp6[site][mu];
						}
					}
					break;
				}
				else r_norm = norm_next;

			}
			
		}*/

		//Now we try to solve the dirac equation
		/*r = source;
		AlgebraUtils::setToZero(solution);
		int conjugateSpaceDimension = 7;
		std::complex<real_t> rho[conjugateSpaceDimension];

		for (int j = 0; j < 535; ++j) {
			//Inner conjugate gradient solver
			{
				multigrid_vector_t solution_hat, r_hat, p_hat, tmp_hat;
				multiGridProjector->apply(solution_hat,r);
				//source_hat = solution_hat;
				multiGridOperator->multiply(tmp_hat,solution_hat);

				//std::cout << toString(source_hat.asVector()) << std::endl;

#pragma omp parallel for
				for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
					r_hat[m] = solution_hat[m] - tmp_hat[m];
					p_hat[m] = r_hat[m];
				}

				long_real_t norm = 0.;
#pragma omp parallel for reduction(+:norm)
				for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
					norm += real(conj(r_hat[m])*r_hat[m]);
				}
				reduceAllSum(norm);
				long_real_t norm_next = norm;

				for (unsigned int innerStep = 0; innerStep < 50; ++innerStep) {
					multiGridOperator->multiply(tmp_hat,p_hat);
					norm = norm_next;
					long_real_t gammaRe = 0;
					long_real_t gammaIm = 0;
#pragma omp parallel for reduction(+:gammaRe,gammaIm)
					for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
						std::complex<real_t> result = conj(p_hat[m])*tmp_hat[m];
						gammaRe += real(result);
						gammaIm += imag(result);
					}
					reduceAllSum(gammaRe);
					reduceAllSum(gammaIm);
					std::complex<real_t> alpha = static_cast<real_t>(norm)/std::complex<real_t>(gammaRe,gammaIm);


#pragma omp parallel for
					for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
						solution_hat[m] = solution_hat[m] + alpha * p_hat[m];
						r_hat[m] = r_hat[m] - alpha * tmp_hat[m];
					}

					norm_next = 0.;
#pragma omp parallel for reduction(+:norm_next)
					for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
						norm_next += real(conj(r_hat[m])*r_hat[m]);
					}
					//Check for convergence
					if (norm_next < 0.000000000001) break;
					//std::cout << "Inner norm at step " << innerStep << ": " << norm_next << std::endl; 
					
					real_t beta = static_cast<real_t>(norm_next/norm);

#pragma omp parallel for
					for (int m = 0; m < multigrid_vector_t::Layout::size; ++m) {
						p_hat[m] = r_hat[m] + beta * p_hat[m];
					}
				}

				//std::cout << toString(solution_hat.asVector()) << std::endl;

				multiGridProjector->apply(w[(j % 15)],solution_hat);
			}

			//w[(j % conjugateSpaceDimension)] = r;

			d[(j % conjugateSpaceDimension)] = w[(j % conjugateSpaceDimension)];

			int js = ((j % conjugateSpaceDimension) == 0) ? 0 : 1;
			for (int k = js; k < (j % conjugateSpaceDimension); ++k) {
				std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(w[(j % conjugateSpaceDimension)],ad[k]))/rho[k];
#pragma omp parallel for
				for (int site = 0; site < guess.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						d[(j % conjugateSpaceDimension)][site][mu] = d[(j % conjugateSpaceDimension)][site][mu] - alpha*d[k][site][mu];
					}
				}
			}

			squareDiracWilsonOperator->multiply(ad[(j % conjugateSpaceDimension)],d[(j % conjugateSpaceDimension)]);
			rho[(j % conjugateSpaceDimension)] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[(j % conjugateSpaceDimension)],ad[(j % conjugateSpaceDimension)]));

			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[(j % conjugateSpaceDimension)],r))/rho[(j % conjugateSpaceDimension)];

#pragma omp parallel for
			for (int site = 0; site < guess.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = solution[site][mu] + proj*d[(j % conjugateSpaceDimension)][site][mu];
					r[site][mu] = r[site][mu] - proj*ad[(j % conjugateSpaceDimension)][site][mu];
				}
			}
			//std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;
			if (AlgebraUtils::squaredNorm(r) < 0.00000000001) {
				stepsMax = j;
				break;
			}
			else std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;

		}*/

		

		multigrid_vector_t tmp1, tmp2;
		for (int i = 0; i < multigrid_vector_t::Layout::size; ++i) tmp2[i] = std::complex<real_t>(0.2-i*i,0.4+i);
		
		int numberTests;
		try {
			numberTests = environment.configurations.get<unsigned int>("number_multiplication_test_speed");
		} catch (NotFoundOption& e) {
			numberTests = 300;
		}

		//struct timespec start, finish;
		//double elapsed;
		clock_gettime(CLOCK_REALTIME, &start);
		for (int i = 0; i < numberTests; ++i) {
			multiGridOperator->multiply(tmp1,tmp2);
		}
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
#ifdef TEST_PAPI_SPEED
		if((retval=PAPI_flops( &real_time, &proc_time, &flpins, &mflops)) < PAPI_OK) test_fail(__FILE__, __LINE__, "PAPI_flops", retval);
		if (isOutputProcess()) std::cout << "MFLOPS for MultiGridOperator: " << mflops << " MFLOPS. " << std::endl;
#endif
		if (isOutputProcess()) std::cout << "Timing for MultiGridOperator: " << (elapsed*1000)/numberTests << " ms."<< std::endl;
		

		delete multiGridOperator;
	}
#endif
	
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
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Gamma5 test on DiracWilsonOperator: " << g5test << std::endl;
		ImprovedDiracWilsonOperator* disastro2 = new ImprovedDiracWilsonOperator();
		disastro2->setLattice(environment.getFermionLattice());
		disastro2->setKappa(0.15);
		disastro2->setCSW(1.);
		disastro2->multiply(test3, test1);
		disastro2->setGamma5(false);
		disastro2->multiply(test2, test1);
		g5test = norm(AlgebraUtils::dot(test3,test3) - AlgebraUtils::gamma5dot(test3,test2));
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Gamma5 test on ImprovedDiracWilsonOperator: " << g5test << std::endl;
	}
	
	/*{
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.000000001);
		biConjugateGradient->setMaximumSteps(300);
		
		BlockDiracWilsonOperator* blackBlockDiracWilsonOperator = new BlockDiracWilsonOperator(Black);
		blackBlockDiracWilsonOperator->setLattice(environment.getFermionLattice());
		blackBlockDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		blackBlockDiracWilsonOperator->setGamma5(false);
		
		BlockDiracWilsonOperator* redBlockDiracWilsonOperator = new BlockDiracWilsonOperator(Red);
		redBlockDiracWilsonOperator->setLattice(environment.getFermionLattice());
		redBlockDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		redBlockDiracWilsonOperator->setGamma5(false);
		
		DiracWilsonOperator* diracWilsonOperator = new DiracWilsonOperator();
		diracWilsonOperator->setLattice(environment.getFermionLattice());
		diracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		diracWilsonOperator->setGamma5(false);

		SquareDiracWilsonOperator* squareDiracWilsonOperator = new SquareDiracWilsonOperator();
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());
		squareDiracWilsonOperator->setKappa(environment.configurations.get<double>("kappa"));
		
		ComplementBlockDiracOperator* K = new ComplementBlockDiracOperator(diracWilsonOperator, redBlockDiracWilsonOperator, blackBlockDiracWilsonOperator);
		
		SAPPreconditioner *preconditioner = new SAPPreconditioner(diracWilsonOperator,K);
		preconditioner->setSteps(9);
		preconditioner->setPrecision(0.0001);
		
		reduced_dirac_vector_t source, guess, tmp1, tmp2, tmp3, tmp4, tmp5, Kb, w[15], d[15], ad[15], r, solution;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);
		
		//preconditioner->setSteps(17);
		//preconditioner->setPrecision(0.0000000001);
		struct timespec start, finish;
		double elapsed;
		clock_gettime(CLOCK_REALTIME, &start);

		r = source;
		AlgebraUtils::setToZero(solution);
		std::complex<real_t> rho[15];

		for (int j = 0; j < 15; ++j) {
			if (j == 0 || j == 1) {
				preconditioner->setSteps(17);
				preconditioner->setPrecision(0.0000000001);

#pragma omp parallel for
				for (int site = 0; site < guess.completesize; ++site) {
					tmp2[site][0] = r[site][0];
					tmp2[site][1] = r[site][1];
					tmp2[site][2] = -r[site][2];
					tmp2[site][3] = -r[site][3];
				}
				preconditioner->multiply(guess,tmp2);
#pragma omp parallel for
				for (int site = 0; site < guess.completesize; ++site) {
					guess[site][2] = -guess[site][2];
					guess[site][3] = -guess[site][3];
				}
				preconditioner->multiply(w[j],guess);
			} 
			else {
				w[j] = r;
			}

			d[j] = w[j];

			int js = (j == 0) ? 0 : 1;
			for (int k = js; k < j; ++k) {
				std::complex<real_t> alpha = static_cast< std::complex<real_t> >(AlgebraUtils::dot(w[j],ad[k]))/rho[k];
#pragma omp parallel for
				for (int site = 0; site < guess.completesize; ++site) {
					for (unsigned int mu = 0; mu < 4; ++mu) {
						d[j][site][mu] = d[j][site][mu] - alpha*d[k][site][mu];
					}
				}
			}

			squareDiracWilsonOperator->multiply(ad[j],d[j]);
			rho[j] = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[j],ad[j]));

			std::complex<real_t> proj = static_cast< std::complex<real_t> >(AlgebraUtils::dot(d[j],r))/rho[j];

#pragma omp parallel for
			for (int site = 0; site < guess.completesize; ++site) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					solution[site][mu] = solution[site][mu] + proj*d[j][site][mu];
					r[site][mu] = r[site][mu] - proj*ad[j][site][mu];
				}
			}
			//std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;
			if (AlgebraUtils::squaredNorm(r) < 0.00000000001) {
				stepsMax = j;
				break;
			}
			else std::cout << "Residual norm at step " << j << ":" << AlgebraUtils::squaredNorm(r) << std::endl;

		}

#pragma omp parallel for
		for (int site = 0; site < guess.completesize; ++site) {
			tmp2[site][0] = source[site][0];
			tmp2[site][1] = source[site][1];
			tmp2[site][2] = -source[site][2];
			tmp2[site][3] = -source[site][3];
		}
		preconditioner->multiply(tmp3,tmp2);
#pragma omp parallel for
		for (int site = 0; site < guess.completesize; ++site) {
			tmp3[site][2] = -tmp3[site][2];
			tmp3[site][3] = -tmp3[site][3];
		}
		preconditioner->multiply(tmp4,tmp3);

		biConjugateGradient->solve(squareDiracWilsonOperator,source,solution,&tmp4);
		
		clock_gettime(CLOCK_REALTIME, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Timing for MSAP preconditioning: " << (elapsed) << " s."<< std::endl;
		int stepsMax = biConjugateGradient->getLastSteps();
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Preconditioned number of steps: " << stepsMax << std::endl;

		diracWilsonOperator->setGamma5(true);
		diracWilsonOperator->multiply(tmp2,solution);
		diracWilsonOperator->multiply(tmp3,tmp2);
		std::cout << "Convergence of: " << AlgebraUtils::differenceNorm(tmp3,source) << std::endl;
		
		//redBlockDiracWilsonOperator->multiply(tmp1,source);
		
		//biConjugateGradient->solve(redBlockDiracWilsonOperator, source, tmp1);
	}*/

	/*{
		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.000000001);
		biConjugateGradient->setMaximumSteps(8000);
		reduced_dirac_vector_t source, tmp1, tmp2, tmp3, tmp4, tmp5;
		//First we take a random vector for tests
		AlgebraUtils::generateRandomVector(source);
		AlgebraUtils::normalize(source);
		
		DiracOperator* squareDiracWilsonOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		squareDiracWilsonOperator->setLattice(environment.getFermionLattice());

		//Time and inversion with the standard BiConjugateGradient
		{
			struct timespec start, finish;
			double elapsed;
			clock_gettime(CLOCK_REALTIME, &start);

			biConjugateGradient->solve(squareDiracWilsonOperator, source, tmp1);

			clock_gettime(CLOCK_REALTIME, &finish);
			elapsed = (finish.tv_sec - start.tv_sec);
			elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Timing for standard BiConjugateGradient: " << (elapsed) << " s."<< std::endl;
			int stepsMax = biConjugateGradient->getLastSteps();
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Undeflated number of steps: " << stepsMax << std::endl;
		}

		DeflationInverter* deflationInverter = new DeflationInverter(1);
		deflationInverter->setPrecision(0.0000000001);
		deflationInverter->setBasisDimension(43);
		deflationInverter->setBlockSize(4);
		deflationInverter->generateBasis(squareDiracWilsonOperator);

		//Time and inversion with the deflated BiConjugateGradient
		{
			struct timespec start, finish;
			double elapsed;
			clock_gettime(CLOCK_REALTIME, &start);

			//deflationInverter->solve(squareDiracWilsonOperator, source, tmp2);

			clock_gettime(CLOCK_REALTIME, &finish);
			elapsed = (finish.tv_sec - start.tv_sec);
			elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Timing for deflated BiConjugateGradient: " << (elapsed) << " s."<< std::endl;
			int stepsMax = biConjugateGradient->getLastSteps();
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Deflated number of steps: " << stepsMax << std::endl;
		}
		
		double normError1 = AlgebraUtils::differenceNorm(tmp1,tmp2);
		squareDiracWilsonOperator->multiply(tmp3, tmp2);
		double normError2 = AlgebraUtils::differenceNorm(tmp3,source);
		squareDiracWilsonOperator->multiply(tmp4, tmp1);
		double normError3 = AlgebraUtils::differenceNorm(tmp4,source);
		if (isOutputProcess()) std::cout << "Ridiamo di gusto: " << normError1 << " " << normError2 << " " << " " << normError3 << " " << deflationInverter->getLastSteps() << std::endl;
		if (isOutputProcess()) std::cout << "Zum beispiel: " << std::endl << tmp2[5][3] - tmp1[5][3] << std::endl;
		delete biConjugateGradient;
		delete squareDiracWilsonOperator;
		delete deflationInverter;
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
	/*{
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
	}*/

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
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on DiracWilsonOperator: " << htest << std::endl;
		AlgebraUtils::generateRandomVector(test1);
		AlgebraUtils::generateRandomVector(test2);
		ImprovedDiracWilsonOperator* disastro2 = new ImprovedDiracWilsonOperator();
		disastro2->setLattice(environment.getFermionLattice());
		disastro2->setKappa(0.1);
		disastro2->setCSW(1.);
		disastro2->multiply(test3,test1);
		disastro2->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on ImprovedDiracWilsonOperator: " << htest << std::endl;
		BlockDiracWilsonOperator* disastro3 = new BlockDiracWilsonOperator();
		disastro3->setLattice(environment.getFermionLattice());
		disastro3->setKappa(0.1);
		disastro3->multiply(test3,test1);
		disastro3->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on BlockDiracWilsonOperator: " << htest << std::endl;
		SquareBlockDiracWilsonOperator* disastro4 = new SquareBlockDiracWilsonOperator();
		disastro4->setLattice(environment.getFermionLattice());
		disastro4->setKappa(0.1);
		disastro4->multiply(test3,test1);
		disastro4->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on SquareBlockDiracWilsonOperator: " << htest << std::endl;
		SquareComplementBlockDiracWilsonOperator* disastro5 = new SquareComplementBlockDiracWilsonOperator();
		disastro5->setLattice(environment.getFermionLattice());
		disastro5->setKappa(0.1);
		disastro5->multiply(test3,test1);
		disastro5->multiply(test4,test2);
		htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
		if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on SquareComplementBlockDiracWilsonOperator: " << htest << std::endl;
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

	environment.gaugeLinkConfiguration.updateHalo();
	environment.synchronize();
}

} /* namespace Update */
