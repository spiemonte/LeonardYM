#include "LandauGaugeFixing.h"
#include "ToString.h"

namespace Update {

LandauGaugeFixing::LandauGaugeFixing() : LatticeSweep(), GaugeFixing() { }

LandauGaugeFixing::LandauGaugeFixing(const LandauGaugeFixing& toCopy) : LatticeSweep(toCopy), GaugeFixing(toCopy) { }

LandauGaugeFixing::~LandauGaugeFixing() { }

void LandauGaugeFixing::generateOverrelaxationTransformation(extended_matrix_lattice_t& gauge_transformation, const extended_gauge_lattice_t& lattice, int d) {
	typedef extended_matrix_lattice_t LT;
	typedef extended_matrix_lattice_t::Layout Layout;

#pragma omp parallel for
	for (int site = 0; site < gauge_transformation.localsize; ++site) {
		if ((Layout::globalIndexX(site) + Layout::globalIndexY(site) + Layout::globalIndexZ(site) + Layout::globalIndexT(site)) % 2 == (d % 2)) {
			GaugeGroup staple;
			set_to_zero(staple);
			for (unsigned int mu = 0; mu < 4; ++mu) {
				staple += lattice[site][mu];
				staple += htrans(lattice[LT::sdn(site,mu)][mu]);
			}
#if NUMCOLORS > 2
			set_to_identity(gauge_transformation[site]);
			for (unsigned int k = 0; k < numberColors-1; ++k) {
				for (int l = k+1; l < numberColors; ++l) {
					GaugeGroup Q;
					set_to_identity(Q);
					Q(k,k) = staple.at(k,k);
					Q(l,l) = staple.at(l,l);
					Q(k,l) = staple.at(k,l);
					Q(l,k) = staple.at(l,k);
					Eigen::ComplexEigenSolver<GaugeGroup> es(htrans(Q)*Q);
					GaugeGroup inverse_square_root;
					inverse_square_root.zeros();
					for (int i = 0; i < numberColors; ++i) {
						inverse_square_root.at(i,i) = 1./sqrt(es.eigenvalues()[i]);
					}
					GaugeGroup m = es.eigenvectors();
					GaugeGroup omega = Q * m * inverse_square_root * htrans(m);
					gauge_transformation[site] = htrans(omega)*gauge_transformation[site];
					staple = htrans(omega)*staple;
				}
			}

			/*GaugeGroup Q, R;
			qr(Q,R,staple);//Compute the QR decomposition
			//Now set the det(Q) to 1:
			std::complex<real_t> determinant = det(Q);
			// by multiply the first row with the complex conjugate of det:
			for (int i = 0; i < numberColors; ++i) {
				for (int j = 0; j < numberColors; ++j) {
					Q(i,j) /= pow(determinant,1./numberColors);
				}
			}
			if (real(trace(htrans(Q)*staple)) > real(trace(staple))) {
				gauge_transformation[site] = htrans(Q);
			}
			else {
				set_to_identity(gauge_transformation[site]);
			}
			if (site == 0) std::cout << "S: "<< toString(staple) << std::endl;
			if (site == 0) std::cout << "S: "<< toString(gauge_transformation[site]) << std::endl;
			if (site == 0) exit(13);*/
/*
#ifdef EIGEN
			Eigen::ComplexEigenSolver<GaugeGroup> es(htrans(staple)*staple);
			GaugeGroup inverse_square_root;
			inverse_square_root.zeros();
			for (int i = 0; i < numberColors; ++i) {
				inverse_square_root.at(i,i) = 1./sqrt(es.eigenvalues()[i]);
			}
			GaugeGroup m = es.eigenvectors();
			gauge_transformation[site] = staple * m * inverse_square_root * htrans(m);
			if (real(trace(gauge_transformation[site]*staple)) < real(trace(staple))) {
				set_to_identity(gauge_transformation[site]);
			}
#endif
*/

#endif
#if NUMCOLORS == 2
			//if (site == 0) std::cout << "S: "<< toString(staple) << std::endl;
			std::complex<real_t> determinant = det(staple);
			for (int i = 0; i < numberColors; ++i) {
				for (int j = 0; j < numberColors; ++j) {
					staple(i,j) /= sqrt(determinant);
				}
			}
			gauge_transformation[site] = htrans(staple);
			//if (site == 0) std::cout << "G: "<< toString(gauge_transformation[site]) << std::endl;
			//if (site == 0) exit(13);
#endif
		}
		else {
			set_to_identity(gauge_transformation[site]);
		}
	}

	gauge_transformation.updateHalo();
}

void LandauGaugeFixing::execute(environment_t& environment) {
	typedef extended_matrix_lattice_t LT;
	typedef extended_matrix_lattice_t::Layout Layout;

	real_t epsilon1 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon1");
	real_t epsilon2 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon2");
	real_t epsilon3 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon3");

	real_t beta1 = environment.configurations.get<real_t>("LandauGaugeFixing::beta1");
	real_t beta2 = environment.configurations.get<real_t>("LandauGaugeFixing::beta2");
	real_t beta3 = environment.configurations.get<real_t>("LandauGaugeFixing::beta3");

	unsigned int steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::steps");

	real_t acceptance1 = 0., acceptance2 = 0., acceptance3 = 0., acceptancept = 0.;

	extended_gauge_lattice_t gaugeFixed1 = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t gaugeFixed2 = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t gaugeFixed3 = environment.gaugeLinkConfiguration;
	extended_gauge_lattice_t maximum     = environment.gaugeLinkConfiguration;

	long_real_t functionalValue1 = functional(gaugeFixed1);
	long_real_t functionalValue2 = functionalValue1;
	long_real_t functionalValue3 = functionalValue2;
	long_real_t maximalFunctionalValue = functionalValue1;

	extended_gauge_lattice_t tmp;

	extended_matrix_lattice_t gauge_transformation;

	for (unsigned int i = 0; i < steps; ++i) {
		//The first
		tmp = gaugeFixed1;
		this->generateRandomGaugeTransformation(gauge_transformation, epsilon1);
		this->transform(tmp, gauge_transformation);
		long_real_t newFunctional = functional(tmp);
		
		if (this->metropolis(beta1*(newFunctional - functionalValue1))) {
			functionalValue1 = newFunctional;
			gaugeFixed1 = tmp;
			acceptance1 += 1.;
		}
		if (newFunctional > maximalFunctionalValue) {
			maximalFunctionalValue = newFunctional;
			maximum = tmp;
		}
				

		//The second
		tmp = gaugeFixed2;
		this->generateRandomGaugeTransformation(gauge_transformation, epsilon2);
		this->transform(tmp, gauge_transformation);
		newFunctional = functional(tmp);
		
		if (this->metropolis(beta2*(newFunctional - functionalValue2))) {
			functionalValue2 = newFunctional;
			gaugeFixed2 = tmp;
			acceptance2 += 1.;
		}
		if (newFunctional > maximalFunctionalValue) {
			maximalFunctionalValue = newFunctional;
			maximum = tmp;
		}
		

		//The third
		tmp = gaugeFixed3;
		this->generateRandomGaugeTransformation(gauge_transformation, epsilon3);
		this->transform(tmp, gauge_transformation);
		newFunctional = functional(tmp);
		
		if (this->metropolis(beta3*(newFunctional - functionalValue3))) {
			functionalValue3 = newFunctional;
			gaugeFixed3 = tmp;
			acceptance3 += 1.;
		}
		if (newFunctional > maximalFunctionalValue) {
			maximalFunctionalValue = newFunctional;
			maximum = tmp;
		}
		

		//Parallel tempering
		if (this->metropolis((beta2-beta3)*(functionalValue2 - functionalValue3))) {
			long_real_t ftmp = functionalValue2;
			functionalValue2 = functionalValue3;
			functionalValue3 = ftmp;

			tmp = gaugeFixed2;
			gaugeFixed2 = gaugeFixed3;
			gaugeFixed3 = tmp;

			acceptancept += 1.;
		}

		if (this->metropolis((beta1-beta2)*(functionalValue1 - functionalValue2))) {
			long_real_t ftmp = functionalValue1;
			functionalValue1 = functionalValue2;
			functionalValue2 = ftmp;

			tmp = gaugeFixed1;
			gaugeFixed1 = gaugeFixed2;
			gaugeFixed2 = tmp;

			acceptancept += 1.;
		}
	}

	if (isOutputProcess()) {
		std::cout << "LandauGaugeFixing::MC Acceptances (global phase): " << acceptance1/steps << " " << acceptance2/steps << " " << acceptance3/steps << std::endl;
		std::cout << "LandauGaugeFixing::PT Acceptances (global phase): " << acceptancept/(2.*steps) << std::endl;
	}

	acceptance1 = 0., acceptance2 = 0., acceptance3 = 0., acceptancept = 0.;

	unsigned int local_steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::local_steps");

	long_real_t convergence = 0;

	for (unsigned int i = 0; i < local_steps; ++i) {
		if ((i) % static_cast<int>(local_steps/(environment.configurations.get<unsigned int>("LandauGaugeFixing::output_steps"))) == 0 && isOutputProcess()) {
			if (isOutputProcess()) std::cout << "LandauGaugeFixing::Maximal functional at step " << i  << ": " << maximalFunctionalValue << std::endl;
			if (i > 0 && isOutputProcess()) std::cout << "LandauGaugeFixing::   convergence: " << convergence << std::endl;
			if (fabs(convergence) < environment.configurations.get<real_t>("LandauGaugeFixing::precision") && i > 0) break;
			convergence = 0;
		}
		
		tmp = maximum;
		this->generateOverrelaxationTransformation(gauge_transformation, tmp, i);
		this->transform(tmp, gauge_transformation);
		long_real_t newFunctional = functional(tmp);
		
		if (newFunctional - maximalFunctionalValue > 0.) {
			convergence += newFunctional - maximalFunctionalValue;
			maximalFunctionalValue = newFunctional;
			maximum = tmp;
		}

		/*tmp = maximum;
		this->generateRandomGaugeTransformation(gauge_transformation, epsilon1);

		for (int site = 0; site < tmp.localsize; ++site) {
			GaugeGroup oldS[4], oldD[4];
			real_t resultOld = 0.;
			for (unsigned int mu = 0; mu < 4; ++mu) {
				resultOld += real(trace(tmp[site][mu]));
				oldS[mu] = tmp[site][mu];
				resultOld += real(trace(tmp[LT::sdn(site,mu)][mu]));
				oldD[mu] = tmp[LT::sdn(site,mu)][mu];
			}
	
			for (unsigned int mu = 0; mu < 4; ++mu) {
				tmp[site][mu] = htrans(gauge_transformation[site])*tmp[site][mu];
				tmp[LT::sdn(site,mu)][mu] = tmp[LT::sdn(site,mu)][mu]*gauge_transformation[site];
			}
	
			real_t resultNew = 0.;
			for (unsigned int mu = 0; mu < 4; ++mu) {
				resultNew += real(trace(tmp[site][mu]));
				resultNew += real(trace(tmp[LT::sdn(site,mu)][mu]));
			}

			if (resultNew < resultOld) {
				for (unsigned int mu = 0; mu < 4; ++mu) {
					tmp[site][mu] = oldS[mu];
					tmp[LT::sdn(site,mu)][mu] = oldD[mu];
				}
			}
			
		}

		newFunctional = functional(tmp);
		if (newFunctional - maximalFunctionalValue > 0.) {
			std::cout << "Giusto per: " << newFunctional - maximalFunctionalValue << std::endl;
			convergence += newFunctional - maximalFunctionalValue;
			maximalFunctionalValue = newFunctional;
			maximum = tmp;
			acceptance1 += 1.;
		}*/
		
	}

	environment.gaugeLinkConfiguration = maximum;
	environment.synchronize();

	if (isOutputProcess()) std::cout << "LandauGaugeFixing::Final deviation from the Landau gauge: " << this->deviation(environment.gaugeLinkConfiguration) << std::endl;
}

long_real_t LandauGaugeFixing::functional(const extended_gauge_lattice_t& lattice) {
	long_real_t result = 0.;
	
#pragma omp parallel for reduction(+:result)
	for (int site = 0; site < lattice.localsize; ++site) {
		for (unsigned int mu = 0; mu < 4; ++mu) {
			result += real(trace(lattice[site][mu]));
		}
	}

	reduceAllSum(result);
	return result;
}

long_real_t LandauGaugeFixing::deviation(const extended_gauge_lattice_t& lattice) const {
	typedef extended_matrix_lattice_t LT;
	typedef extended_matrix_lattice_t::Layout Layout;

	extended_gauge_lattice_t Afield;
	this->getLieAlgebraField(Afield, lattice);

	long_real_t convergence = 0.;

#pragma omp parallel for reduction(+:convergence)
	for (int site = 0; site < Afield.localsize; ++site) {
		GaugeGroup dmu;
		set_to_zero(dmu);
		for (unsigned int mu = 0; mu < 4; ++mu) {
			dmu += (Afield[site][mu]-Afield[LT::sdn(site,mu)][mu]);
		}
		convergence += real(trace(dmu*htrans(dmu)));
	}

	reduceAllSum(convergence);
	return convergence/static_cast<real_t>(Layout::globalVolume);
}

void LandauGaugeFixing::registerParameters(po::options_description& desc) {
	desc.add_options()
		("LandauGaugeFixing::epsilon1", po::value<real_t>()->default_value(0.000001), "set the epsilon for the random transformation used by the first PT run")
		("LandauGaugeFixing::epsilon2", po::value<real_t>()->default_value(0.00000001), "set the epsilon for the random transformation used by the second PT run")
		("LandauGaugeFixing::epsilon3", po::value<real_t>()->default_value(0.00000000001), "set the epsilon for the random transformation used by the third PT run")

		("LandauGaugeFixing::beta1", po::value<real_t>()->default_value(5.), "set the beta for the random transformation used by the first PT run")
		("LandauGaugeFixing::beta2", po::value<real_t>()->default_value(1.), "set the beta for the random transformation used by the second PT run")
		("LandauGaugeFixing::beta3", po::value<real_t>()->default_value(0.2), "set the beta for the random transformation used by the third PT run")

		("LandauGaugeFixing::steps", po::value<unsigned int>()->default_value(500), "set the number of MC trials")
		("LandauGaugeFixing::local_steps", po::value<unsigned int>()->default_value(3000), "set the number of local over-relaxation trials")
		("LandauGaugeFixing::precision", po::value<real_t>()->default_value(0.000000000001), "Set the covergence precision")
		("LandauGaugeFixing::output_steps", po::value<unsigned int>()->default_value(100), "Set the step to monitor the output")
		;
}

}

