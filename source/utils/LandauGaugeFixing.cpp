#include "LandauGaugeFixing.h"

namespace Update {

LandauGaugeFixing::LandauGaugeFixing() : LatticeSweep(), GaugeFixing() { }

LandauGaugeFixing::LandauGaugeFixing(const LandauGaugeFixing& toCopy) : LatticeSweep(toCopy), GaugeFixing(toCopy) { }

LandauGaugeFixing::~LandauGaugeFixing() { }

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

	for (unsigned int i = 0; i < 5*steps; ++i) {
		if ((i+1) % steps == 0 && isOutputProcess()) std::cout << "LandauGaugeFixing::Maximal functional at step " << i + 1 << ": " << maximalFunctionalValue << std::endl;
		//The first
		tmp = gaugeFixed1;
		this->generateRandomGaugeTransformation(gauge_transformation, epsilon1);
#ifdef ENABLE_MPI
		for (int processor = 0; processor < Layout::numberProcessors; ++processor) {
			if (processor == Layout::this_processor) {
#endif
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
#ifdef ENABLE_MPI			
			}
			else {
				for (int site = tmp.sharedsize; site < tmp.localsize; ++site) {
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
			}
			tmp.updateHalo();
		}
#endif

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
	}

	environment.gaugeLinkConfiguration = maximum;
	environment.synchronize();
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

void LandauGaugeFixing::registerParameters(po::options_description& desc) {
	desc.add_options()
		("LandauGaugeFixing::epsilon1", po::value<real_t>()->default_value(0.000001), "set the epsilon for the random transformation used by the first PT run")
		("LandauGaugeFixing::epsilon2", po::value<real_t>()->default_value(0.00000001), "set the epsilon for the random transformation used by the second PT run")
		("LandauGaugeFixing::epsilon3", po::value<real_t>()->default_value(0.00000000001), "set the epsilon for the random transformation used by the third PT run")

		("LandauGaugeFixing::beta1", po::value<real_t>()->default_value(5.), "set the beta for the random transformation used by the first PT run")
		("LandauGaugeFixing::beta2", po::value<real_t>()->default_value(1.), "set the beta for the random transformation used by the second PT run")
		("LandauGaugeFixing::beta3", po::value<real_t>()->default_value(0.2), "set the beta for the random transformation used by the third PT run")

		("LandauGaugeFixing::steps", po::value<unsigned int>()->default_value(500), "set the number of MC trials")
		;
}

}

