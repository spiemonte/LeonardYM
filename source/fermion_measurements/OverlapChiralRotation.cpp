#include "OverlapChiralRotation.h"
#include "io/GlobalOutput.h"
#include "algebra_utils/AlgebraUtils.h"
#include "utils/StoutSmearing.h"
#include "utils/Gamma.h"
#include "inverters/PreconditionedBiCGStab.h"
#include "dirac_operators/Propagator.h"

namespace Update {

OverlapChiralRotation::OverlapChiralRotation() : diracOperator(0), inverter(0) { }

OverlapChiralRotation::OverlapChiralRotation(const OverlapChiralRotation& toCopy) : LatticeSweep(toCopy), StochasticEstimator(), diracOperator(0), inverter(0) { }

OverlapChiralRotation::~OverlapChiralRotation() {
	if (diracOperator) delete diracOperator;
	if (inverter) delete inverter;
}

void OverlapChiralRotation::chiral_rotation_psi(DiracOperator* dirac, const reduced_dirac_vector_t& input, reduced_dirac_vector_t& output, unsigned int steps, const std::complex<real_t>& theta) const {
	output = input;
	reduced_dirac_vector_t tmp1, tmp2;
	
	for (unsigned int i = 0; i < steps; ++i) {
		tmp1 = output;
		dirac->multiply(tmp2, output);

#pragma omp parallel for
		for (int site = 0; site < output.completesize; ++site) {
			for (unsigned int mu = 0; mu < 2; ++mu) {
				output[site][mu] +=  (theta/static_cast<real_t>(steps))*tmp1[site][mu] - (theta/static_cast<real_t>(steps))*tmp2[site][mu];
			}
			for (unsigned int mu = 2; mu < 4; ++mu) {
				output[site][mu] +=  -(theta/static_cast<real_t>(steps))*tmp1[site][mu] + (theta/static_cast<real_t>(steps))*tmp2[site][mu];
			}
		}
	}
}

void OverlapChiralRotation::chiral_rotation_psibar(DiracOperator* dirac, const reduced_dirac_vector_t& input, reduced_dirac_vector_t& output, unsigned int steps, const std::complex<real_t>& theta) const {
	output = input;
	reduced_dirac_vector_t tmp1, tmp2;
	
	for (unsigned int i = 0; i < steps; ++i) {
		tmp1 = output;
		AlgebraUtils::gamma5(tmp1);
		dirac->multiply(tmp2, tmp1);

#pragma omp parallel for
		for (int site = 0; site < output.completesize; ++site) {
			for (unsigned int mu = 0; mu < 4; ++mu) {
				output[site][mu] +=  (theta/static_cast<real_t>(steps))*tmp1[site][mu] - (theta/static_cast<real_t>(steps))*tmp2[site][mu];
			}
		}
	}
}

void OverlapChiralRotation::execute(environment_t& environment) {
	typedef extended_dirac_vector_t::Layout Layout;//TODO: only vector operations?

	//unsigned int max_step = environment.configurations.get<unsigned int>("number_stochastic_estimators");
	if (diracOperator == 0) {
		diracOperator = DiracOperator::getInstance(environment.configurations.get<std::string>("dirac_operator"), 1, environment.configurations);
	}

	extended_fermion_lattice_t lattice;

	try {
		unsigned int numberLevelSmearing = environment.configurations.get<unsigned int>("OverlapChiralRotation::levels_stout_smearing");
		double smearingRho = environment.configurations.get<double>("OverlapChiralRotation::rho_stout_smearing");
		StoutSmearing stoutSmearing;
#ifdef ADJOINT
		extended_gauge_lattice_t smearedConfiguration;
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, smearedConfiguration, numberLevelSmearing, smearingRho);
		ConvertLattice<extended_fermion_lattice_t,extended_gauge_lattice_t>::convert(lattice, smearedConfiguration);//TODO
#endif
#ifndef ADJOINT
		stoutSmearing.spatialSmearing(environment.gaugeLinkConfiguration, lattice, numberLevelSmearing, smearingRho);
#endif
		environment.setFermionBc(lattice);
	} catch (NotFoundOption& ex) {
		lattice =  environment.getFermionLattice();
		if (isOutputProcess()) std::cout << "OverlapChiralRotation::No smearing options found, proceeding without!" << std::endl;
	}

	try {
		unsigned int t = environment.configurations.get<unsigned int>("OverlapChiralRotation::t_source_origin");
		extended_fermion_lattice_t swaplinkconfig;
		typedef extended_fermion_lattice_t LT;
		if (t != 0) {
			for (unsigned int n = 0; n < t; ++n) {
				//We do a swap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) swaplinkconfig[site][mu] = lattice[site][mu];
				}
				swaplinkconfig.updateHalo();
				//We wrap
				for(int site = 0; site < (lattice.localsize); ++site){
					for (unsigned int mu = 0; mu < 4; ++mu) lattice[site][mu] = swaplinkconfig[LT::sup(site,3)][mu];
				}
				lattice.updateHalo();
			}
		}
	} catch (NotFoundOption& ex) {

	}

	diracOperator->setLattice(lattice);
	diracOperator->setGamma5(false);
	
	if (inverter == 0) inverter = new BiConjugateGradient();
	inverter->setPrecision(environment.configurations.get<double>("OverlapChiralRotation::inverter_precision"));
	inverter->setMaximumSteps(environment.configurations.get<unsigned int>("OverlapChiralRotation::inverter_max_steps"));

	std::vector< std::complex<long_real_t> > unrotated_chiral_condensate;
	std::vector< std::complex<long_real_t> > rotated_chiral_condensate;

	unsigned int numberStochasticEstimators = environment.configurations.get<unsigned int>("OverlapChiralRotation::number_stochastic_estimators");
	unsigned int expProductsteps = environment.configurations.get<unsigned int>("OverlapChiralRotation::exponential_product_number_steps");

	int inversionSteps = 0;
	for (unsigned int i = 0; i < numberStochasticEstimators; ++i) {
		this->generateRandomNoise(randomNoise);
		
		inverter->solve(diracOperator, randomNoise, inverse);
		Propagator::constructPropagator(diracOperator, inverse, tmp);

		std::complex<long_real_t> condensate = AlgebraUtils::dot(randomNoise, tmp);

		//Rotate the fermion fields
		unrotated_chiral_condensate.push_back(condensate);

		if (isOutputProcess()) std::cout << "OverlapChiralRotation::Unrotated condensate for source " << i << " " << condensate << std::endl;

		this->chiral_rotation_psi(diracOperator, inverse, psi_rotated, expProductsteps, std::complex<real_t>(0.,PI/2.));

		Propagator::constructPropagator(diracOperator, psi_rotated, tmp);

		this->chiral_rotation_psibar(diracOperator, tmp, barpsi_rotated, expProductsteps, std::complex<real_t>(0.,PI/2.));

		condensate = AlgebraUtils::dot(randomNoise, barpsi_rotated);

		rotated_chiral_condensate.push_back(condensate);

		if (isOutputProcess()) std::cout << "OverlapChiralRotation::Rotated condensate for source " << i << " " << condensate << std::endl;
	}

	if (isOutputProcess()) std::cout << "OverlapChiralRotation::residual mass computed with " << inversionSteps << " inversion steps" << std::endl;
	
	if (environment.measurement && isOutputProcess()) {
		GlobalOutput* output = GlobalOutput::getInstance();

		output->push("unrotated_chiral_condensate");

		std::complex<long_real_t> mean_unrotated_chiral_condensate = this->mean(unrotated_chiral_condensate);

		output->write("unrotated_chiral_condensate", real(mean_unrotated_chiral_condensate));
		output->write("unrotated_chiral_condensate", imag(mean_unrotated_chiral_condensate));

		std::cout << "OverlapChiralRotation::unrotated_chiral_condensate " << mean_unrotated_chiral_condensate << " +/- " << this->standardDeviation(unrotated_chiral_condensate)/static_cast<long double>(sqrt(numberStochasticEstimators)) << std::endl;
		
		output->pop("unrotated_chiral_condensate");

		output->push("rotated_chiral_condensate");

		std::complex<long_real_t> mean_rotated_chiral_condensate = this->mean(rotated_chiral_condensate);

		output->write("rotated_chiral_condensate", real(mean_rotated_chiral_condensate));
		output->write("rotated_chiral_condensate", imag(mean_rotated_chiral_condensate));

		std::cout << "OverlapChiralRotation::rotated_chiral_condensate " << mean_rotated_chiral_condensate << " +/- " << this->standardDeviation(rotated_chiral_condensate)/static_cast<long double>(sqrt(numberStochasticEstimators)) << std::endl;
	}
}

void OverlapChiralRotation::registerParameters(po::options_description& desc) {
	desc.add_options()
		("OverlapChiralRotation::number_stochastic_estimators", po::value<unsigned int>()->default_value(20), "The number of stochastic estimators to be used")
		("OverlapChiralRotation::exponential_product_number_steps", po::value<unsigned int>()->default_value(200), "The number of steps to be used to approximate the chiral rotation")
		("OverlapChiralRotation::inverter_precision", po::value<real_t>()->default_value(0.0000000001), "set the inverter precision")
		("OverlapChiralRotation::inverter_max_steps", po::value<unsigned int>()->default_value(10000), "maximum number of inverter steps")
		("OverlapChiralRotation::rho_stout_smearing", po::value<real_t>(), "set the stout smearing parameter")
		("OverlapChiralRotation::levels_stout_smearing", po::value<unsigned int>(), "levels of stout smearing")
		;
}

} /* namespace Update */
