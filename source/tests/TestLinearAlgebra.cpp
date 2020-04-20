#include "TestLinearAlgebra.h"
#include "inverters/BiConjugateGradient.h"
#include "inverters/GMRESR.h"
#include "algebra_utils/AlgebraUtils.h"
#include "dirac_operators/SquareDiracWilsonOperator.h"
#include "dirac_operators/SquareImprovedDiracWilsonOperator.h"
#include "dirac_operators/DiracWilsonOperator.h"
#include "dirac_operators/ImprovedDiracWilsonOperator.h"
#include "dirac_operators/OverlapOperator.h"
#include "dirac_operators/ExactOverlapOperator.h"
#include "dirac_operators/BasicDiracWilsonOperator.h"
#include "dirac_operators/SquareBlockDiracWilsonOperator.h"
#include "dirac_operators/ComplementBlockDiracOperator.h"
#include "dirac_operators/SquareComplementBlockDiracWilsonOperator.h"
#include "dirac_operators/SquareComplementBlockDiracOperator.h"
#include "dirac_operators/SquareTwistedDiracOperator.h"
#include "dirac_operators/TwistedDiracOperator.h"
#include "dirac_operators/SAPPreconditioner.h"
#include "multigrid/MultiGridOperator.h"
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
		try {
			OverlapOperator* disastro6 = new OverlapOperator();
			disastro6->setLattice(environment.getFermionLattice());
			disastro6->setKappa(0.12);
			disastro6->setMass(0.1);

			std::vector< std::complex<real_t> > coeff = environment.configurations.get< std::vector< std::complex<real_t> > >("OverlapOperator::squareRootApproximation");
			Polynomial squareRootApproximation;
			squareRootApproximation.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			squareRootApproximation.setRoots(coeff);

			disastro6->setSquareRootApproximation(squareRootApproximation);

			disastro6->multiply(test3,test1);
			disastro6->multiply(test4,test2);
			htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on OverlapOperator: " << htest << std::endl;
			delete disastro6;

			ExactOverlapOperator* disastro7 = new ExactOverlapOperator();
			disastro7->setLattice(environment.getFermionLattice());
			disastro7->setKappa(0.2);
			disastro7->setMass(0.1);

			disastro7->setSquareRootApproximation(squareRootApproximation);

			disastro7->multiply(test3,test1);
			disastro7->multiply(test4,test2);
			htest = norm(AlgebraUtils::dot(test2,test3) - AlgebraUtils::dot(test4,test1));
			if (isOutputProcess()) std::cout << "TestLinearAlgebra::Hermitian test on ExactOverlapOperator: " << htest << std::endl;
			delete disastro7;
		}
		catch (Update::NotFoundOption& e) {
			std::cout << "TestLinearAlgebra::Proceeding without Hermitian test on Overlap operator" << std::endl;
		}
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
