/*
 * NFlavorQCDUpdater.cpp
 *
 *  Created on: May 2, 2012
 *      Author: spiem_01
 */

#include "NFlavorBlockUpdater.h"
#include "NFlavorQCDAction.h"
#include "AlgebraUtils.h"
#include "MultishiftSolver.h"
#include "RationalApproximation.h"
#include "BiConjugateGradient.h"
#include "DiracEigenSolver.h"
#include "GaugeAction.h"
#include "Integrate.h"
#include "ToString.h"
#include "./dirac_operators/BlockDiracOperator.h"
#include "./dirac_operators/SquareComplementBlockDiracWilsonOperator.h"
#include "./dirac_operators/SquareDiracWilsonOperator.h"
#include <iomanip>

//#define DEBUGFORCE
//#define REVERSIBILITY_CHECK

#ifdef DEBUGFORCE
#include "TestForce.h"
#endif

namespace Update {

NFlavorBlockUpdater::NFlavorBlockUpdater() : NFlavorQCDUpdater() { }

NFlavorBlockUpdater::NFlavorBlockUpdater(const NFlavorBlockUpdater& toCopy) : NFlavorQCDUpdater(toCopy) { }

NFlavorBlockUpdater::~NFlavorBlockUpdater() { }

void NFlavorBlockUpdater::initializeCorrectionStepApproximations(const environment_t& environment) {
	int breakupLevel = environment.configurations.get< unsigned int >("correction_step_breakup_level");
	
	if (polynomialApproximationsInverse.empty()) {
		for (int i = 1; i <= breakupLevel; ++i) {
			std::vector< std::complex<real_t> > coeff = environment.configurations.get< std::vector< std::complex<real_t> > >(std::string("correction_approximation_inverse_")+toString(i));
			Polynomial approx;
			approx.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			approx.setRoots(coeff);
			polynomialApproximationsInverse.push_back(approx);
		}
	}

	if (polynomialApproximationsDirect.empty()) {
		for (int i = 1; i <= breakupLevel; ++i) {
			std::vector< std::complex<real_t> > coeff = environment.configurations.get< std::vector< std::complex<real_t> > >(std::string("correction_approximation_direct_")+toString(i));
			Polynomial approx;
			approx.setScaling(coeff.front());
			coeff.erase(coeff.begin(),coeff.begin()+1);
			approx.setRoots(coeff);
			polynomialApproximationsDirect.push_back(approx);
		}
	}
}

void NFlavorBlockUpdater::checkCorrectionStepApproximations(const environment_t& environment) const {
	double testerDirect = 1., testerInverse = 1.;
	double epsilon = 0.0001;
	int breakupLevel = environment.configurations.get< unsigned int >("correction_step_breakup_level");
	for (int i = 0; i < breakupLevel; ++i) {
		testerInverse *= polynomialApproximationsInverse[i].evaluate(2.).real();
		testerDirect *= polynomialApproximationsDirect[i].evaluate(2.).real();
	}
	std::cout << "Vediamo da dove partire: " << testerInverse << " " << testerDirect << std::endl;
	double numberFermions = -2*log(testerInverse)/log(2.);//We are using the square of the dirac operator, hence we have a factor 2
	if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::The theory has " <<  numberFermions << " nf." << std::endl;
#ifdef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 0.5) < epsilon) std::cout << "NFlavorQCDUpdater::The theory seems SUSY" << std::endl;
#endif
#ifndef ADJOINT
	if (isOutputProcess() && fabs(numberFermions - 1) < epsilon) std::cout << "NFlavorBlockUpdater::The theory seems 1 FlavorQCD" << std::endl;
	if (isOutputProcess() && fabs(numberFermions - 2) < epsilon) std::cout << "NFlavorBlockUpdater::The theory seems 2 FlavorQCD" << std::endl;
#endif
	if (isOutputProcess() && fabs(testerInverse*testerDirect - 1.) > epsilon) {
		std::cout << "NFlavorBlockUpdater::Warning, large mismatch between force and metropolis approximations: " << fabs(testerInverse*testerDirect - 1.) << std::endl;
	}
}

void NFlavorBlockUpdater::execute(environment_t& environment) {	
	//First we initialize the approximations
	this->initializeApproximations(environment);

	//Then we initialize the approximations for the correction step
	this->initializeCorrectionStepApproximations(environment);

	//We check the theory that is simulated
	if (environment.iteration == 0 && environment.sweep == 0) {
		this->checkTheory(environment);
		this->checkCorrectionStepApproximations(environment);
	}

	//Initialize the momenta
	this->randomMomenta(momenta);
	//Copy the environment
	environmentNew = environment;

	//Take the Block Dirac Operator
	if (diracOperator == 0)  diracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	diracOperator->setLattice(environment.getFermionLattice());

	//Take the Block Square Dirac Operator
	if (squareDiracOperator == 0) squareDiracOperator = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
	squareDiracOperator->setLattice(environment.getFermionLattice());

	//Take the gauge action
	if (gaugeAction == 0) gaugeAction = GaugeAction::getInstance(environment.configurations.get<std::string>("name_action"),environment.configurations.get<double>("beta"));

	//Initialize the pseudofermion fields
	std::vector<extended_dirac_vector_t>::iterator i;
	std::vector<RationalApproximation>::iterator j = rationalApproximationsHeatBath.begin();
	long_real_t oldPseudoFermionEnergy = 0.;
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		this->generateGaussianDiracVector(tmp_pseudofermion);
		oldPseudoFermionEnergy += AlgebraUtils::squaredNorm(tmp_pseudofermion);
		//Heat bath by multiply the tmp_pseudofermion to the polynomial of the dirac operator
		j->evaluate(squareDiracOperator, *i, tmp_pseudofermion);//TODO
		++j;
	}

	if (environment.sweep == 0 && environment.iteration == 0 && environment.measurement) {
		//Now we test the correctness of rational/heatbath
		long_real_t test = 0.;
		//Now we use the better approximation for the metropolis step
		std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
		for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
			//Now we evaluate it with the rational approximation of the inverse
			rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
			test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
			++rational;
		}
		if (isOutputProcess()) std::cout << "NFlavoBlockUpdater::Consistency check for the metropolis: " << test - oldPseudoFermionEnergy << std::endl;
		//Now we use the approximation for the force step
		rational = rationalApproximationsForce.begin();
		test = 0.;
		for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
			//Now we evaluate it with the rational approximation of the inverse
			rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
			test += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
			++rational;
		}
		if (isOutputProcess()) std::cout << "NFlavoBlockUpdater::Consistency check for the force: " << test - oldPseudoFermionEnergy << std::endl;
	}

	//Get the initial energy of momenta
	long_real_t oldMomentaEnergy = this->momentaEnergy(momenta);
	//Get the initial energy of the lattice
	long_real_t oldLatticeEnergy = gaugeAction->energy(environment);

	//Set the dirac operator to the new environment
	squareDiracOperator->setLattice(environmentNew.getFermionLattice());
	diracOperator->setLattice(environmentNew.getFermionLattice());

	//Take the fermion action
	if (fermionAction == 0) fermionAction = new NFlavorFermionAction(squareDiracOperator, diracOperator, rationalApproximationsForce);
	//TODO dirac operator correctly set?
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		//Add the pseudofermions
		fermionAction->addPseudoFermion(&(*i));
	}
	//Set the precision of the inverter
	fermionAction->setForcePrecision(environment.configurations.get<double>("force_inverter_precision"));

	//Take the global action
	if (nFlavorQCDAction == 0) nFlavorQCDAction = new NFlavorAction(gaugeAction, fermionAction);

	//The t-length of a single integration step
	real_t t_length = environment.configurations.get<double>("hmc_t_length");
	std::vector<Force*> forces;
	//The numbers of integration steps
	std::vector<unsigned int> numbers_steps = environment.configurations.get< std::vector<unsigned int> >("number_hmc_steps");
	if (numbers_steps.size() == 1) {
		forces.push_back(nFlavorQCDAction);
	} else if (numbers_steps.size() == 2) {
		forces.push_back(fermionAction);
		forces.push_back(gaugeAction);
	}
	else {
		if (isOutputProcess()) std::cout << "NFlavorBlockUpdater::Warning, NFlavor does not support more than two time integrations!" << std::endl;
		numbers_steps.resize(1);
		forces.push_back(nFlavorQCDAction);
	}

#ifdef DEBUGFORCE
	//Test of the force
	TestForce testForce;
	extended_gauge_lattice_t tmp;
	//Internal calculations needed
	nFlavorQCDAction->updateForce(tmp, environment);
	testForce.genericTestForce(environment, nFlavorQCDAction, tmp[5][2], 5, 2);
#endif

	Integrate* integrate = Integrate::getInstance(environment.configurations.get<std::string>("name_integrator"));

	//Integrate numerically the equation of motion
	integrate->integrate(environmentNew, momenta, forces, numbers_steps, t_length);
#ifdef REVERSIBILITY_CHECK
	integrate->integrate(environmentNew, momenta, forces, numbers_steps, -t_length);
#endif

	//Get the final energy of momenta
	long_real_t newMomentaEnergy = this->momentaEnergy(momenta);
	//Get the final energy of the lattice
	long_real_t newLatticeEnergy = gaugeAction->energy(environmentNew);
	//Get the final energy of the pseudofermions
	long_real_t newPseudoFermionEnergy = 0.;
	diracOperator->setLattice(environmentNew.getFermionLattice());//TODO TODO TODO
	squareDiracOperator->setLattice(environmentNew.getFermionLattice());
	//Now we use the better approximation for the metropolis step
	std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		//Now we evaluate it with the rational approximation of the inverse
		rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
		newPseudoFermionEnergy += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
		++rational;
	}

	//Global Metropolis Step
	bool metropolis = this->metropolis(oldMomentaEnergy + oldLatticeEnergy + oldPseudoFermionEnergy, newMomentaEnergy + newLatticeEnergy + newPseudoFermionEnergy);
	if (metropolis) {
		this->stochasticCorrectionStep(environment, environmentNew);
		//Now we do the correction step
		/*DiracOperator* squareDiracOperatorCorrection = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);
		DiracOperator* squareBlockDiracOperatorCorrection = BlockDiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 2, environment.configurations);

		squareDiracOperatorCorrection->setLattice(environment.getFermionLattice());
		squareBlockDiracOperatorCorrection->setLattice(environment.getFermionLattice());

		long_real_t oldCorrectionEnergy = 0., newCorrectionEnergy = 0.;

		//First we generate and heatbath the gaussian noise vector
		reduced_dirac_vector_t gaussianVectors[polynomialApproximationsInverse.size()];
		reduced_dirac_vector_t noiseVector, tmp;

		std::vector<Polynomial>::iterator directApproximation = polynomialApproximationsDirect.begin();
		std::vector<Polynomial>::iterator inverseApproximation = polynomialApproximationsInverse.begin();
		
		for (unsigned int i = 0; i < polynomialApproximationsInverse.size(); ++i) {
			this->generateGaussianDiracVector(noiseVector);
			oldCorrectionEnergy += AlgebraUtils::squaredNorm(noiseVector);
			//Heat bath by multiply the noiseVector to the polynomial of the dirac operator
			//directApproximation->evaluate(squareDiracOperatorCorrection, tmp, noiseVector);
			//inverseApproximation->evaluate(squareBlockDiracOperatorCorrection, gaussianVectors[i], tmp);
			directApproximation->evaluate(squareDiracOperatorCorrection, gaussianVectors[i], noiseVector);
			++directApproximation;
			++inverseApproximation;
		}

		//Now we evaluate the new correction factor energy
		//squareDiracOperatorCorrection->setLattice(environmentNew.getFermionLattice());
		//squareBlockDiracOperatorCorrection->setLattice(environmentNew.getFermionLattice());

		directApproximation = polynomialApproximationsDirect.begin();
		inverseApproximation = polynomialApproximationsInverse.begin();
		for (unsigned int i = 0; i < polynomialApproximationsInverse.size(); ++i) {
			//The energy now is computed in the reverse order!
			//directApproximation->evaluate(squareBlockDiracOperatorCorrection, tmp, gaussianVectors[i]);
			//Use noiseVector as a temporary storage
			//inverseApproximation->evaluate(squareDiracOperatorCorrection, noiseVector, tmp);
			inverseApproximation->evaluate(squareDiracOperatorCorrection, noiseVector, gaussianVectors[i]);
			++directApproximation;
			++inverseApproximation;
			newCorrectionEnergy += AlgebraUtils::squaredNorm(noiseVector);
		}

		std::cout << "Giusto per: " << oldCorrectionEnergy << " " << newCorrectionEnergy << std::endl;
		
		bool correction = this->metropolis(oldCorrectionEnergy, newCorrectionEnergy);
		if (correction) {
			environment = environmentNew;
	}

	delete squareDiracOperatorCorrection;
	delete squareBlockDiracOperatorCorrection;*/

		/*std::cout << "Giusto per: " << newPseudoFermionEnergy - oldPseudoFermionEnergy << std::endl;

		BiConjugateGradient* biConjugateGradient = new BiConjugateGradient();
		biConjugateGradient->setPrecision(0.0001);
		reduced_dirac_vector_t random, temp, result;
		AlgebraUtils::generateRandomGaussianVector(temp);*/

		/*ComplementBlockDiracWilsonOperator* complementBlockDiracWilsonOperatorOld = new ComplementBlockDiracWilsonOperator();
		complementBlockDiracWilsonOperatorOld->setKappa(environment.configurations.get<real_t>("kappa"));
		complementBlockDiracWilsonOperatorOld->setLattice(environment.getFermionLattice());
		complementBlockDiracWilsonOperatorOld->setPrecision(0.000001);

		ComplementBlockDiracWilsonOperator* complementBlockDiracWilsonOperatorNew = new ComplementBlockDiracWilsonOperator();
		complementBlockDiracWilsonOperatorNew->setKappa(environment.configurations.get<real_t>("kappa"));
		complementBlockDiracWilsonOperatorNew->setLattice(environmentNew.getFermionLattice());
		complementBlockDiracWilsonOperatorNew->setPrecision(0.000001);*/

		/*DiracWilsonOperator* diracWilsonOperatorOld = new DiracWilsonOperator();
		diracWilsonOperatorOld->setKappa(environment.configurations.get<real_t>("kappa"));
		diracWilsonOperatorOld->setLattice(environment.getFermionLattice());

		DiracWilsonOperator* diracWilsonOperatorNew = new DiracWilsonOperator();
		diracWilsonOperatorNew->setKappa(environment.configurations.get<real_t>("kappa"));
		diracWilsonOperatorNew->setLattice(environmentNew.getFermionLattice());

		long_real_t mean = 0.;
		long_real_t meansq = 0.;

		for (int i = 0; i < 10; ++i) {
			AlgebraUtils::generateRandomGaussianVector(random);
			
			biConjugateGradient->solve(diracWilsonOperatorOld, random, temp);
			diracWilsonOperatorNew->multiply(result, temp);
			

			real_t partial = real(AlgebraUtils::dot(random, random) - AlgebraUtils::dot(result,result));
			mean += partial;
			meansq += partial*partial;
			
			std::cout << "Partial result: " << partial << std::endl;
		}

		std::cout << "NFlavorBlockUpdate::Stochastic mean: " << mean/10 << std::endl;
		std::cout << "NFlavorBlockUpdate::Stochastic standard deviation: " << sqrt((meansq/10.) - (mean/10.)*(mean/10.))/sqrt(10.) << std::endl;*/
		/*long_real_t meanNew = 0.;

		complementBlockDiracWilsonOperator->setLattice(environmentNew.getFermionLattice());

		for (int i = 0; i < 30; ++i) {
			AlgebraUtils::generateRandomGaussianVector(temp);
			complementBlockDiracWilsonOperator->multiply(result, temp);

			meanNew += real(AlgebraUtils::dot(temp, temp) - AlgebraUtils::dot(result,result));
		}
		std::cout << "Media stocastica: " << meanNew/30 << std::endl;*/

		/*std::cout << "NFlavorBlockUpdate::Determinant ratio: " << sqrt(exp(-mean/10.)) << std::endl;	
		real_t determinantRatio = sqrt(exp(-mean/10));
		bool correction = this->metropolis(determinantRatio);*/
		
		/*if (correction) {
			environment = environmentNew;
		}*/
	}

	//Determinant test
	/*{
		reduced_dirac_vector_t temp, result;
		AlgebraUtils::generateRandomGaussianVector(temp);

		ComplementBlockDiracWilsonOperator* complementBlockDiracWilsonOperator = new ComplementBlockDiracWilsonOperator();
		complementBlockDiracWilsonOperator->setKappa(environment.configurations.get<real_t>("kappa"));
		complementBlockDiracWilsonOperator->setLattice(environment.getFermionLattice());

		long_real_t meanOld = 0.;

		for (int i = 0; i < 100; ++i) {
			AlgebraUtils::generateRandomGaussianVector(temp);
			complementBlockDiracWilsonOperator->multiply(result, temp);

			meanOld += real(AlgebraUtils::dot(temp, temp) - AlgebraUtils::dot(result,result));
		}

		std::cout << "Media stocastica: " << meanOld/100 << std::endl;
		long_real_t meanNew = 0.;

		complementBlockDiracWilsonOperator->setLattice(environmentNew.getFermionLattice());

		for (int i = 0; i < 100; ++i) {
			AlgebraUtils::generateRandomGaussianVector(temp);
			complementBlockDiracWilsonOperator->multiply(result, temp);

			meanNew += real(AlgebraUtils::dot(temp, temp) - AlgebraUtils::dot(result,result));
		}
		std::cout << "Media stocastica: " << meanNew/100 << std::endl;

		std::cout << "Determinant ratio: " << sqrt(exp(-(meanNew/100) + (meanOld/100))) << std::endl;
	}*/


	delete integrate;
}

void NFlavorBlockUpdater::stochasticCorrectionStep(environment_t& environment, const environment_t& environmentNew) {
	DiracOperator* squareDiracOperator = new SquareComplementBlockDiracWilsonOperator();//
	squareDiracOperator->setKappa(environment.configurations.get<real_t>("kappa"));
	squareDiracOperator->setLattice(environment.getFermionLattice());
	
	//Initialize the pseudofermion fields
	std::vector<extended_dirac_vector_t>::iterator i;
	std::vector<RationalApproximation>::iterator j = rationalApproximationsHeatBath.begin();
	long_real_t oldPseudoFermionEnergy = 0.;
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		this->generateGaussianDiracVector(tmp_pseudofermion);
		oldPseudoFermionEnergy += AlgebraUtils::squaredNorm(tmp_pseudofermion);
		//Heat bath by multiply the tmp_pseudofermion to the polynomial of the dirac operator
		j->evaluate(squareDiracOperator, *i, tmp_pseudofermion);//TODO
		++j;
	}

	squareDiracOperator->setLattice(environmentNew.getFermionLattice());
	
	long_real_t newPseudoFermionEnergy = 0.;
	//Now we use the better approximation for the metropolis step
	std::vector<RationalApproximation>::iterator rational = rationalApproximationsMetropolis.begin();
	for (i = pseudofermions.begin(); i != pseudofermions.end(); ++i) {
		//Now we evaluate it with the rational approximation of the inverse
		rational->evaluate(squareDiracOperator,tmp_pseudofermion,*i);
		newPseudoFermionEnergy += real(AlgebraUtils::dot(*i,tmp_pseudofermion));
		++rational;
	}

	bool metropolis = this->metropolis(oldPseudoFermionEnergy, newPseudoFermionEnergy);
	if (metropolis) {
		environment = environmentNew;
	}
}

} /* namespace Update */
