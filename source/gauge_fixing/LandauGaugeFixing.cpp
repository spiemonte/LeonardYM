#include "LandauGaugeFixing.h"
#include "utils/ToString.h"

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
					if (fabs(trace(htrans(omega)*omega) - (double)numberColors) > 1e-10) {
						continue;
					}
					gauge_transformation[site] = htrans(omega)*gauge_transformation[site];
					staple = htrans(omega)*staple;
				}
			}

#endif
#if NUMCOLORS == 2
			std::complex<real_t> determinant = det(staple);
			for (int i = 0; i < numberColors; ++i) {
				for (int j = 0; j < numberColors; ++j) {
					staple(i,j) /= sqrt(determinant);
				}
			}
			gauge_transformation[site] = htrans(staple);
#endif
		}
		else {
			set_to_identity(gauge_transformation[site]);
		}
	}

	gauge_transformation.updateHalo();
}

void LandauGaugeFixing::execute(environment_t& environment) {
	real_t epsilon1 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon1");
	real_t epsilon2 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon2");
	real_t epsilon3 = environment.configurations.get<real_t>("LandauGaugeFixing::epsilon3");

	real_t beta1 = environment.configurations.get<real_t>("LandauGaugeFixing::beta1");
	real_t beta2 = environment.configurations.get<real_t>("LandauGaugeFixing::beta2");
	real_t beta3 = environment.configurations.get<real_t>("LandauGaugeFixing::beta3");

	unsigned int steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::steps");
	unsigned int local_steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::local_steps");

	real_t precision = environment.configurations.get<real_t>("LandauGaugeFixing::precision");

	unsigned int final_steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::final_steps");
	unsigned int number_copies = environment.configurations.get<unsigned int>("LandauGaugeFixing::number_copies");

	unsigned int output_steps = environment.configurations.get<unsigned int>("LandauGaugeFixing::output_steps");

	std::vector<long_real_t> maximal_values(number_copies);
	std::vector<extended_gauge_lattice_t> maximals(number_copies);
	for (unsigned int i = 0; i < number_copies; ++i) {
		maximals[i] = environment.gaugeLinkConfiguration;
		maximal_values[i] = this->gaugeFixing(maximals[i], epsilon1, beta1, epsilon2, beta2, epsilon3, beta3, steps, local_steps, precision, output_steps);
		if (isOutputProcess()) std::cout << "LandauGaugeFixing::Maximal functional for copy " << i << ": " << maximal_values[i] << std::endl;
	}

	long_real_t maximum = 0.;
	unsigned int i_maximum = 0;
	for (unsigned int i = 0; i < number_copies; ++i) {
		if (maximal_values[i] > maximum) {
			maximum = maximal_values[i];
			i_maximum = i;
		}
	}

	this->gaugeFixing(maximals[i_maximum], epsilon1, beta1, epsilon2, beta2, epsilon3, beta3, steps, final_steps, precision, output_steps);
	environment.gaugeLinkConfiguration = maximals[i_maximum];
	environment.synchronize();
	
	long_real_t deviation = this->deviation(environment.gaugeLinkConfiguration);

	if (isOutputProcess()) std::cout << "LandauGaugeFixing::Final deviation from the Landau gauge: " << deviation << std::endl;
}

long_real_t LandauGaugeFixing::gaugeFixing(extended_gauge_lattice_t& lattice, const real_t& epsilon1, const real_t& beta1, const real_t& epsilon2, const real_t& beta2, const real_t& epsilon3, const real_t& beta3, unsigned int steps, unsigned int local_steps, const real_t& precision, unsigned int output_steps) {
	typedef extended_matrix_lattice_t LT;
	typedef extended_matrix_lattice_t::Layout Layout;

	real_t acceptance1 = 0., acceptance2 = 0., acceptance3 = 0., acceptancept = 0.;

	extended_gauge_lattice_t gaugeFixed1 = lattice;
	extended_gauge_lattice_t gaugeFixed2 = lattice;
	extended_gauge_lattice_t gaugeFixed3 = lattice;
	extended_gauge_lattice_t maximum     = lattice;

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

	long_real_t convergence = 0;

	for (unsigned int i = 0; i < local_steps; ++i) {
		if ((i) % static_cast<int>(local_steps/output_steps) == 0) {
			if (isOutputProcess()) std::cout << "LandauGaugeFixing::Maximal functional at step " << i  << ": " << maximalFunctionalValue << std::endl;
			if (i > 0 && isOutputProcess()) std::cout << "LandauGaugeFixing::   convergence: " << convergence << std::endl;
			if (fabs(convergence) < precision && i > 0) break;
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
	}

	lattice = maximum;
	return maximalFunctionalValue;
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

void LandauGaugeFixing::registerParameters(std::map<std::string, Option>& desc) {
	desc["LandauGaugeFixing::epsilon1"] = Option("LandauGaugeFixing::epsilon1", 1e-5, "set the epsilon for the random transformation used by the first PT run");
	desc["LandauGaugeFixing::epsilon2"] = Option("LandauGaugeFixing::epsilon2", 1e-7, "set the epsilon for the random transformation used by the second PT run");
	desc["LandauGaugeFixing::epsilon3"] = Option("LandauGaugeFixing::epsilon3", 1e-9, "set the epsilon for the random transformation used by the third PT run");

	desc["LandauGaugeFixing::beta1"] = Option("LandauGaugeFixing::beta1", 5.0, "set the beta for the random transformation used by the first PT run");
	desc["LandauGaugeFixing::beta2"] = Option("LandauGaugeFixing::beta2", 1.0, "set the beta for the random transformation used by the second PT run");
	desc["LandauGaugeFixing::beta3"] = Option("LandauGaugeFixing::beta3", 0.2, "set the beta for the random transformation used by the third PT run");

	desc["LandauGaugeFixing::steps"] = Option("LandauGaugeFixing::steps", 500, "set the number of MC trials");
	desc["LandauGaugeFixing::local_steps"] = Option("LandauGaugeFixing::local_steps", 3000, "set the number of local over-relaxation trials");
	desc["LandauGaugeFixing::final_steps"] = Option("LandauGaugeFixing::final_steps", 3000, "set the number of final local over-relaxation trials");
	desc["LandauGaugeFixing::precision"] = Option("LandauGaugeFixing::precision", 1e-12, "Set the covergence precision");
	desc["LandauGaugeFixing::output_steps"] = Option("LandauGaugeFixing::output_steps", 100, "Set the step to monitor the output");
	desc["LandauGaugeFixing::number_copies"] = Option("LandauGaugeFixing::number_copies", 5, "Set the of copies of maxima to use to search the global maximum");
}

}

